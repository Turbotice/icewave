"""
tsync_v3.py -- clock synchronisation laptop <-> phone over UDP (merges tsync.py and tsync_windows.py fro mphonefleet).

PRINCIPLE
  The laptop sends a 4-byte UDP request (value 1) to the phone (port 5000); the phone answers on port 5001 with its
  own clock (8 bytes). With t1 = time before sending, t2 = after sending, t3 = reception, the lag laptop-phone is
  estimated as the median of (t1 - t_phone + t3 - t_phone)/2 over the fastest exchanges (round trip <= median).
  The phone must first be told to listen: GET http://<ip>:8080/udp-sync (done by time_sync_ip).

WHAT CHANGED vs tsync.py / tsync_windows.py
  One file for Linux, macOS and Windows (SO_REUSEADDR on the receiving socket; '<=' median selection, which also
  works on Windows where timer resolution can give identical durations). Removed: per-packet debug prints, the
  matplotlib import, the meaningless socket.socket(SOL_SOCKET, SO_REUSEADDR, 1) line, the hard-coded save folder,
  the unused Windows get_adress() hack. time_sync() now RETURNS its result.

DEPENDENCIES   numpy, stdlib; connect_v3.py, run_gobannos_v3.py, rw_data_v3.py (only for the command line).
USED BY        mobile_v3.py (command 'time')

EXAMPLES
  py tsync_v3.py --phone 36                       # 3 syncs of phone 36, results saved in Phonefleet/Gobannos_Tsync
  py tsync_v3.py --phone 36 --n 200 --iter 5 --folder D:/sync
  import tsync_v3 as ts
  Dt = ts.time_sync_ip("10.153.247.36", n=100, timeout=0.1)   # raw exchanges, or None if unreachable
  ts.get_lag(Dt)   # {'tlag': ..., 'dtmedian': ..., 'tmin': ..., 'tmax': ..., 'tstd': ..., 'n': ..., 't0': ...}
"""
import argparse
import socket
import time
import urllib.request
from pathlib import Path

import numpy as np

PORT = 8080            # Gobannos HTTP port
UDP_PHONE = 5000       # phone listens here
UDP_LOCAL = 5001       # laptop listens here


def time_sync_ip(ip, n=1000, timeout=0.1):
    """Run n UDP exchanges with the phone. Returns {1,2,3: arrays of lags, 'duration': array, 'time': t0} or None."""
    do_nothing = (0).to_bytes(4, "big")
    respond = (1).to_bytes(4, "big")
    stop = (2).to_bytes(4, "big")

    sock_send = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    sock_receive = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    sock_receive.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
    try:
        sock_receive.bind(("", UDP_LOCAL))
        sock_receive.settimeout(timeout)
        sock_send.sendto(do_nothing, (ip, UDP_PHONE))
        try:
            urllib.request.urlopen(f"http://{ip}:{PORT}/udp-sync", timeout=5).read()    # phone starts listening
        except OSError:
            print(f"Unable to connect to {ip}")
            return None

        Dt = {1: [], 2: [], 3: []}
        duration, lost = [], 0
        t0 = time.time()
        for _ in range(n):
            t1 = time.time_ns()
            sock_send.sendto(respond, (ip, UDP_PHONE))
            t2 = time.time_ns()
            try:
                t_phone_bytes = sock_receive.recv(8)
            except OSError:                      # timeout
                lost += 1
                continue
            t3 = time.time_ns()
            t_phone = int.from_bytes(t_phone_bytes, byteorder="big")
            duration.append((t3 - t1) * 1e-9)
            Dt[1].append((t1 - t_phone) * 1e-9)
            Dt[2].append((t2 - t_phone) * 1e-9)
            Dt[3].append((t3 - t_phone) * 1e-9)
        sock_send.sendto(stop, (ip, UDP_PHONE))
        print(f"Duration : {np.round((time.time() - t0) * 1000, 3)} ms; packets lost: {lost}/{n}")
    finally:
        sock_send.close()
        sock_receive.close()
    if not duration:
        print("No answer from the phone.")
        return None
    out = {k: np.asarray(v) for k, v in Dt.items()}
    out["time"] = t0
    out["duration"] = np.asarray(duration)
    return out


def get_lag(Dt):
    """Robust lag summary from time_sync_ip() output (fastest half of the exchanges)."""
    duration = np.asarray(Dt["duration"])
    median_duration = np.median(duration)
    idx = np.where(duration <= median_duration)[0]
    lag = (np.asarray(Dt[1])[idx] + np.asarray(Dt[3])[idx]) / 2
    return {"tlag": np.median(lag), "dtmedian": median_duration, "tmin": np.min(lag), "tmax": np.max(lag),
            "tstd": np.std(lag), "n": len(duration), "t0": Dt["time"]}


def time_sync(phone, n=1000, timeout=0.1, network=None):
    """Sync one phone by Phonefleet ID; returns the raw exchanges (or None)."""
    import connect_v3 as connect
    ip = connect.get_adress(phone, network=network if network is not None else connect.NETWORK)
    return time_sync_ip(ip, n=n, timeout=timeout)


def main():
    import connect_v3 as connect
    import rw_data_v3 as rw
    ap = argparse.ArgumentParser(description="Clock synchronisation laptop <-> phone")
    ap.add_argument("--phone", type=int, default=connect.DEFAULT_PHONELIST[0], help="Phonefleet ID")
    ap.add_argument("--network", type=int, default=connect.NETWORK, help="3rd IP octet")
    ap.add_argument("--n", type=int, default=100, help="exchanges per synchronisation")
    ap.add_argument("--iter", type=int, default=3, help="number of synchronisations")
    ap.add_argument("--folder", default="Phonefleet/Gobannos_Tsync")
    args = ap.parse_args()
    for _ in range(args.iter):
        Dt = time_sync(args.phone, n=args.n, network=args.network)
        if Dt is None:
            continue
        res = {args.phone: get_lag(Dt)}
        print(res)
        rw.writedict_csv(str(Path(args.folder) / f"tsync_{int(round(time.time()))}.txt"), res)


if __name__ == "__main__":
    main()