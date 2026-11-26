"""
run_gobannos_v3.py -- minimal HTTP client for the Gobannos app running on the phone (stdlib only).

The phone serves a small HTTP server (default port 8080). Every function below is one GET request:

    /status              -> 'STOPPED' / 'RUNNING'...        get_status(ip)
    /start?name          -> start the acquisition           individual_start(ip)
    /stop                -> stop the acquisition            individual_stop(ip)
    /list-files          -> '[file1, file2, ...]'           get_file_list(ip)
    /get-file/<name>     -> file content (bytes)            get_file(ip, name)
    /usb-cmd/<cmd>       -> text sent to the Arduino        usb_cmd(ip, cmd)   e.g. 'mo', 'c50.0', 'x'
    /clean-files         -> erase ALL files on the phone    individual_clean(ip)   (asks y/n first)
    /udp-sync            -> enable UDP clock sync           (used by tsync_v3.py directly)

CHANGES vs run_gobannos.py (Stéphane, phonefleet version)
  The v1 file also held a legacy asyncio/aiohttp/phyphox/zip code path:  removed here. time_sync() lives in tsync_v3.py.
  Fixed: usb_cmd() now RAISES when the phone cannot be reached (v1 printed a message and the caller believed the
  command was sent); individual_clean() now asks y/n BEFORE erasing (v1 erased first, then asked);
  get_file_list() no longer returns [''] when the phone has no file; get_status() returns text, not bytes.

DEPENDENCIES   none (stdlib).
USED BY        chipiron_session_v3.py, mobile_v3.py, tsync_v3.py.  The run / oscillate scripts do NOT use it
               (they talk to the phone through chipiron_helpers_v3.HttpLink).

EXAMPLES
    import run_gobannos_v3 as gob
    gob.get_status("10.153.247.36")            # 'STOPPED'
    gob.individual_start("10.153.247.36")
    gob.usb_cmd("10.153.247.36", "mo")          # Arduino -> oscillation mode
    gob.usb_cmd("10.153.247.36", "c60.0")       # 60 % duty
    files = gob.get_file_list("10.153.247.36")
    data = gob.get_file("10.153.247.36", files[0])
"""
import http.client
import urllib.request
from urllib.parse import quote

port = 8080                       # Gobannos HTTP port (same as connect_v3.PORT)
TIMEOUT = 10                      # seconds, for every request

# accessible functions: name -> path
commands = {"start": "/start?name",
            "stop": "/stop",
            "status": "/status",
            "sync": "/kick-sync?threshold",
            "clean": "/clean-files"}


def _open(ip, path, timeout=TIMEOUT):
    """GET http://<ip>:<port><path>; returns the body as bytes."""
    with urllib.request.urlopen(f"http://{ip}:{port}{path}", timeout=timeout) as resp:
        return resp.read()


def get_base_url(ip):
    """'http://<ip>:<port>'."""
    return f"http://{ip}:{port}"


def get_status(ip):
    """Gobannos status as text ('STOPPED', 'RUNNING'...)."""
    status = _open(ip, commands["status"]).decode("utf-8", errors="replace").strip()
    print(f"Status: {status}")
    return status


def get_file_list(ip):
    """List of the files on the phone (the server answers '[a, b, c]')."""
    s = _open(ip, "/list-files").decode("utf-8", errors="replace").strip()
    return [f for f in s[1:-1].split(", ") if f]


def get_file(ip, filename):
    """Content of one file (bytes). If the transfer is cut, the partial content is returned with a warning."""
    try:
        return _open(ip, "/get-file/" + filename, timeout=60)
    except http.client.IncompleteRead as exc:
        print(f"WARNING: file {filename} incomplete ({len(exc.partial)} bytes received)")
        return exc.partial


def usb_cmd(ip, cmd):
    """Send <cmd> to the Arduino through the phone. Raises RuntimeError if the phone does not answer.
    ',' and '.' stay literal (the phone server does not decode '%2C')."""
    try:
        return _open(ip, "/usb-cmd/" + quote(cmd, safe=",."), timeout=5).decode("utf-8", errors="replace").strip()
    except OSError as exc:
        raise RuntimeError(f"not able to send usb instruction {cmd!r} to {ip}: {exc}") from exc


def individual_start(ip):
    """Start the acquisition on the phone."""
    _open(ip, commands["start"])


def individual_stop(ip):
    """Stop the acquisition on the phone."""
    _open(ip, commands["stop"])


def individual_clean(ip):
    """Erase ALL files on the phone, after a y/n confirmation."""
    if input(f"Erase ALL files on the phone {ip}? y/n: ").strip().lower() == "y":
        _open(ip, commands["clean"])
        print("files deleted")
    else:
        print("no delete")