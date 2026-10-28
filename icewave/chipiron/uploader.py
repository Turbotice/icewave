"""Phone-side uploader: scans a folder, uploads finished files, tracks state in SQLite.

Env vars: UPLOAD_TOKEN (required), SERVER_URL, DEVICE_ID, DATA_DIR, STATE_DB, FILE_GLOB
Run:  python uploader.py
"""
import hashlib
import logging
import os
import sqlite3
import time
from pathlib import Path
from urllib.parse import quote

import requests

DATA_DIR = Path(os.environ.get("DATA_DIR", "./recordings"))
STATE_DB = os.environ.get("STATE_DB", "./uploader_state.db")
SERVER_URL = os.environ.get("SERVER_URL", "https://myserver.example.com").rstrip("/")
DEVICE_ID = os.environ.get("DEVICE_ID", "phone-01")
TOKEN = os.environ["UPLOAD_TOKEN"]
FILE_GLOB = os.environ.get("FILE_GLOB", "**/*.csv")
INTERVAL_S = 600          # run every 10 minutes
STABLE_AFTER_S = 60       # a file untouched for 60 s is considered "closed"

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("uploader")


def init_db() -> sqlite3.Connection:
    db = sqlite3.connect(STATE_DB)
    db.execute(
        """CREATE TABLE IF NOT EXISTS files (
               rel_path    TEXT PRIMARY KEY,
               size        INTEGER,
               mtime       REAL,
               sha256      TEXT,
               status      TEXT NOT NULL DEFAULT 'pending',  -- pending | uploaded
               attempts    INTEGER NOT NULL DEFAULT 0,
               last_error  TEXT,
               uploaded_at REAL
           )"""
    )
    db.commit()
    return db


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def scan(db: sqlite3.Connection) -> None:
    """Register new/changed stable files as pending."""
    now = time.time()
    for p in DATA_DIR.glob(FILE_GLOB):
        if not p.is_file():
            continue
        st = p.stat()
        if now - st.st_mtime < STABLE_AFTER_S:
            continue  # still being written
        rel = p.relative_to(DATA_DIR).as_posix()
        row = db.execute("SELECT size, mtime FROM files WHERE rel_path=?", (rel,)).fetchone()
        if row and row == (st.st_size, st.st_mtime):
            continue  # unchanged since last time we looked
        db.execute(
            """INSERT INTO files (rel_path, size, mtime, sha256, status)
               VALUES (?, ?, ?, ?, 'pending')
               ON CONFLICT(rel_path) DO UPDATE SET
                 size=excluded.size, mtime=excluded.mtime,
                 sha256=excluded.sha256, status='pending'""",
            (rel, st.st_size, st.st_mtime, sha256_of(p)),
        )
    db.commit()


def upload_pending(db: sqlite3.Connection) -> None:
    rows = db.execute(
        "SELECT rel_path, sha256 FROM files WHERE status='pending' ORDER BY mtime"
    ).fetchall()
    if not rows:
        log.info("nothing to upload")
        return
    s = requests.Session()
    s.headers["Authorization"] = f"Bearer {TOKEN}"
    for rel, sha in rows:
        url = f"{SERVER_URL}/files/{quote(DEVICE_ID)}/{quote(rel)}"
        try:
            with open(DATA_DIR / rel, "rb") as f:
                r = s.put(url, data=f, headers={"X-SHA256": sha}, timeout=(10, 300))
            r.raise_for_status()
            if r.json().get("sha256") != sha:
                raise RuntimeError("server hash mismatch")
            db.execute(
                "UPDATE files SET status='uploaded', uploaded_at=?, last_error=NULL WHERE rel_path=?",
                (time.time(), rel),
            )
            log.info("uploaded %s", rel)
        except Exception as e:  # network down, file vanished, server error...
            db.execute(
                "UPDATE files SET attempts=attempts+1, last_error=? WHERE rel_path=?",
                (str(e)[:500], rel),
            )
            log.warning("failed %s: %s", rel, e)
            if isinstance(e, requests.ConnectionError):
                break  # offline: stop and retry next cycle
        db.commit()


def main() -> None:
    db = init_db()
    while True:
        started = time.monotonic()
        try:
            scan(db)
            upload_pending(db)
        except Exception:
            log.exception("cycle failed")
        time.sleep(max(1, INTERVAL_S - (time.monotonic() - started)))


if __name__ == "__main__":
    main()
