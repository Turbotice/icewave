"""
mobile_v3.py -- interactive Phonefleet controller (menu).

Controls Gobannos phones over the network: status, start / stop acquisition, list and download files, send USB
commands to the Arduino, measure clock synchronisation.  Use it for one-off actions; chipiron_session_v3.py
chains the same steps automatically for a whole measurement.

DEPENDENCIES   connect_v3.py (IPs), run_gobannos_v3.py (HTTP), tsync_v3.py (clock sync), rw_data_v3.py, numpy
               (tsync_windows.py no longer exists: tsync_v3.py works on every OS).
CONFIGURE      IP_PREFIX / NETWORK / DEFAULT_PHONELIST in connect_v3.py. Run from this folder:  py mobile_v3.py

MENU  (type the command at the '>' prompt)
  network                    change the 3rd IP octet                         > network      then  247
  phones                     change the active phones: 36 or 36:38 (range)   > phones      then  36:38
  status                     status of every active phone                    > status
  time                       clock sync of every active phone (5 rounds)     > time
  start / stop               start / stop the acquisition                    > start
  ls <phone_id>              list files on a phone                           > ls 36
  pull <phone_id> <index>    download one file                               > pull 36 0
  pull <phone_id> <a:b>      download files a..b (inclusive)                 > pull 36 0:4
  usb <command>              send a text command to the Arduino (all active phones; PWM_full_v3.ino):
        usb mo               oscillation mode (rest duty 50 %, bipolar .edc)   usb mc   constant mode (rest 10 %)
        usb c60.0            60 % duty            usb o50,20,0.2   sine 50 +/- 20 % at 0.2 Hz (mode mo only)
        usb x                back to the rest duty                              usb e / usb d   enable / disable motor
  exit
Workflow: status -> start -> experiment -> stop -> ls ID -> pull ID INDEX.
Downloaded files go into OUTPUT_FOLDER as P<phone>_D<original name>.
"""
import os
import sys
import time
from pathlib import Path

import numpy as np

import connect_v3 as connect
import run_gobannos_v3 as gob
import rw_data_v3
import tsync_v3

# ===================== USER CONFIGURATION =====================
OUTPUT_FOLDER = "Phonefleet"  # change to an absolute path if desired
TIME_SYNC_ITERATIONS = 5
# =============================================================


def get_folder():
    """Create and return the output folder with a trailing separator."""
    folder = Path(OUTPUT_FOLDER)
    folder.mkdir(parents=True, exist_ok=True)
    return str(folder) + os.sep


def set_network():
    """Prompt for the third octet of the phone IP address."""
    while True:
        raw = input(f"Network number [{connect.NETWORK}]: ").strip()
        if not raw:
            return connect.NETWORK
        try:
            value = int(raw)
            if 0 <= value <= 255:
                return value
        except ValueError:
            pass
        print("Enter an integer from 0 to 255.")


def set_phonelist():
    """Accept a single phone ID (36) or an inclusive range (36:38)."""
    raw = input(f"Phone ID/range [{':'.join(map(str, connect.DEFAULT_PHONELIST))}]: ").strip()
    if not raw:
        return list(connect.DEFAULT_PHONELIST)
    try:
        parts = [int(part.strip()) for part in raw.split(":")]
        if len(parts) == 1:
            phones = parts
        else:
            phones = []
            for start, end in zip(parts[:-1], parts[1:]):
                step = 1 if end >= start else -1
                phones.extend(range(start, end + step, step))
        if not phones:
            raise ValueError
        return sorted(set(phones))
    except ValueError:
        print("Invalid input. Use e.g. 36, 36:36, or 36:38.")
        return list(connect.DEFAULT_PHONELIST)


def get_filelist(network, phonelist, command, display=True):
    """List remote Gobannos files for a selected phone."""
    parts = command.split()
    if len(parts) < 2:
        print("Usage: ls <phone_id>")
        return []
    try:
        phone = int(parts[1])
    except ValueError:
        print("Phone ID must be an integer.")
        return []
    if phone not in phonelist:
        print(f"Phone {phone} not in active phone list {phonelist}")
        return []
    files = gob.get_file_list(connect.get_adress(phone, network=network))
    if display:
        for index, filename in enumerate(files):
            print(index, filename)
    return files


def save_file(filename, data, phone):
    """Write a downloaded Gobannos file to disk."""
    name = os.path.basename(filename)
    savepath = os.path.join(get_folder(), f"P{phone}_D{name}")
    print(f"Saving: {savepath}")
    lines = data.splitlines()
    print(lines[0] if lines else "(empty file)")
    with open(savepath, "w", encoding="utf-8", newline="") as stream:
        stream.write(data)


def get_file(network, phonelist, command):
    """Download one file or an inclusive range of file indexes."""
    parts = command.split()
    if len(parts) < 3:
        print("Usage: pull <phone_id> <index> or pull <phone_id> <first:last>")
        return
    try:
        phone = int(parts[1])
        selector = parts[2]
        if ":" in selector:
            first, last = map(int, selector.split(":", 1))
            if last < first:
                raise ValueError
            indexes = range(first, last + 1)
        else:
            indexes = [int(selector)]
    except ValueError:
        print("Invalid pull arguments.")
        return
    if phone not in phonelist:
        print(f"Phone {phone} not in active phone list {phonelist}")
        return
    ip = connect.get_adress(phone, network=network)
    files = gob.get_file_list(ip)
    if not files:
        print(f"No files returned by phone {phone} ({ip}).")
        return
    for index in indexes:
        if index < 0 or index >= len(files):
            print(f"Index {index} outside valid range 0..{len(files) - 1}")
            continue
        filename = files[index]
        try:
            raw = gob.get_file(ip, filename)
            data = raw.decode("utf-8") if isinstance(raw, bytes) else str(raw)
            save_file(filename, data, phone)
        except Exception as exc:
            print(f"Could not retrieve {filename}: {exc}")


def check_status(network, phonelist):
    """Query status for all selected phones."""
    if not phonelist:
        print("Phone list is empty; use 'phones' first.")
    for phone in phonelist:
        ip = connect.get_adress(phone, network=network)
        try:
            print(f"Phone {phone} ({ip}): {gob.get_status(ip)}")
        except OSError as exc:
            print(f"Phone {phone} ({ip}) unreachable: {exc}")


def stop_phonelist(network, phonelist):
    """Stop acquisition on all selected phones."""
    for phone in phonelist:
        ip = connect.get_adress(phone, network=network)
        try:
            if gob.get_status(ip) != "STOPPED":
                gob.individual_stop(ip)
        except OSError as exc:
            print(f"Phone {phone} ({ip}): {exc}")


def start_phonelist(network, phonelist):
    """Start acquisition on all selected phones."""
    for phone in phonelist:
        ip = connect.get_adress(phone, network=network)
        try:
            print(f"Phone {phone} ({ip}) before start: {gob.get_status(ip)}")
            gob.individual_start(ip)
        except OSError as exc:
            print(f"Phone {phone} ({ip}): {exc}")


def usb_cmd(network, phonelist, command):
    """Send a USB command to all selected phones."""
    parts = command.split(maxsplit=1)
    if len(parts) < 2:
        print("Usage: usb <command>")
        return
    for phone in phonelist:
        ip = connect.get_adress(phone, network=network)
        try:
            gob.usb_cmd(ip, parts[1])
            print(f"Phone {phone}: sent {parts[1]!r}")
        except RuntimeError as exc:
            print(exc)


def time_sync(network, phonelist, iterations=TIME_SYNC_ITERATIONS):
    """Repeat clock synchronisation and save the results."""
    folder = os.path.join(get_folder(), "Gobannos_Tsync")
    os.makedirs(folder, exist_ok=True)
    for _ in range(iterations):
        results = {}
        for phone in phonelist:
            ip = connect.get_adress(phone, network=network)
            sync_data = tsync_v3.time_sync_ip(ip, n=100, timeout=0.1)
            if sync_data is not None:
                results[phone] = tsync_v3.get_lag(sync_data)
                print(phone, results[phone])
        filename = os.path.join(folder, f"tsync_{int(round(time.time()))}.txt")
        print(f"Saving: {filename}")
        rw_data_v3.writedict_csv(filename, results)


def choose(network, phonelist):
    """Display the command menu and execute one command."""
    print("\nChoose action:")
    for action in ("network", "phones", "status", "time", "start", "stop", "ls <phone_id>",
                   "pull <phone_id> <index|first:last>", "usb <command>", "exit"):
        print(" ", action)
    command = input("> ").strip()
    if command == "network":
        return command, set_network()
    if command == "phones":
        return command, set_phonelist()
    try:
        if command == "status":
            check_status(network, phonelist)
        elif command == "time":
            time_sync(network, phonelist)
        elif command == "start":
            if phonelist:
                start_phonelist(network, phonelist)
            else:
                print("Phone list is empty; use 'phones' first.")
        elif command == "stop":
            stop_phonelist(network, phonelist)
        elif command.startswith("ls"):
            get_filelist(network, phonelist, command)
        elif command.startswith("pull"):
            get_file(network, phonelist, command)
        elif command.startswith("usb"):
            usb_cmd(network, phonelist, command)
        elif command == "exit":
            print("Exiting Phonefleet.")
        elif command:
            print("Command not recognized.")
    except OSError as exc:
        print(f"Network error: {exc}")
    return command, None


def main():
    """Start the interactive controller with the settings of connect_v3.py."""
    network = connect.NETWORK
    phonelist = list(connect.DEFAULT_PHONELIST)
    print("Phonefleet controller")
    print(f"IP scheme: {connect.IP_PREFIX}")
    print(f"Initial network: {network}")
    print(f"Initial phone list: {phonelist}")
    print("Use 'status' first to verify the Gobannos connection.")
    command = ""
    while command != "exit":
        command, output = choose(network, phonelist)
        if command == "phones":
            phonelist = output
        elif command == "network":
            network = output
        print(f"Active phones: {phonelist}")
        print(f"Network: {network}")


if __name__ == "__main__":
    main()