"""Receiving server.

Env vars: UPLOAD_TOKEN (required), STORAGE_DIR
Run:  pip install fastapi uvicorn
      uvicorn server:app --host 0.0.0.0 --port 8000   (put it behind HTTPS: Caddy / nginx)
"""
import hashlib
import hmac
import os
import re
import sqlite3
import time
from pathlib import Path

from fastapi import FastAPI, Header, HTTPException, Request

STORAGE = Path(os.environ.get("STORAGE_DIR", "./storage")).resolve()
TOKEN = os.environ["UPLOAD_TOKEN"]
STORAGE.mkdir(parents=True, exist_ok=True)
INDEX_DB = STORAGE / "index.db"

DEVICE_RE = re.compile(r"^[A-Za-z0-9_-]{1,64}$")
SHA_RE = re.compile(r"^[0-9a-f]{64}$")

app = FastAPI()


def db() -> sqlite3.Connection:
    c = sqlite3.connect(INDEX_DB)
    c.execute(
        """CREATE TABLE IF NOT EXISTS received (
               device_id TEXT, rel_path TEXT, sha256 TEXT, size INTEGER, received_at REAL,
               PRIMARY KEY (device_id, rel_path))"""
    )
    return c


@app.put("/files/{device_id}/{rel_path:path}")
async def put_file(
    device_id: str,
    rel_path: str,
    request: Request,
    x_sha256: str = Header(...),
    authorization: str = Header(""),
):
    if not hmac.compare_digest(authorization, f"Bearer {TOKEN}"):
        raise HTTPException(401, "bad token")
    sha = x_sha256.lower()
    if not DEVICE_RE.match(device_id) or not SHA_RE.match(sha):
        raise HTTPException(400, "bad device id or hash")

    base = (STORAGE / device_id).resolve()
    dest = (base / rel_path).resolve()
    if base not in dest.parents:  # blocks ../ path traversal
        raise HTTPException(400, "bad path")

    conn = db()
    row = conn.execute(
        "SELECT sha256 FROM received WHERE device_id=? AND rel_path=?", (device_id, rel_path)
    ).fetchone()
    if row and row[0] == sha and dest.exists():  # idempotent retry
        return {"status": "already_present", "sha256": sha}

    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_name(dest.name + ".part")
    h, size = hashlib.sha256(), 0
    with open(tmp, "wb") as f:
        async for chunk in request.stream():
            h.update(chunk)
            f.write(chunk)
            size += len(chunk)
        f.flush()
        os.fsync(f.fileno())
    if h.hexdigest() != sha:
        tmp.unlink(missing_ok=True)
        raise HTTPException(422, "hash mismatch, upload corrupted")

    os.replace(tmp, dest)  # atomic: readers never see a partial file
    conn.execute(
        "INSERT OR REPLACE INTO received VALUES (?, ?, ?, ?, ?)",
        (device_id, rel_path, sha, size, time.time()),
    )
    conn.commit()
    conn.close()
    return {"status": "stored", "sha256": sha, "size": size}
