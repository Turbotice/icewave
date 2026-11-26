"""
connect_v3.py -- Phonefleet network configuration: phone ID -> IP address, optional nmap scan.

DEPENDENCIES   stdlib only (nmap program only for scan()/connect()).
USED BY        chipiron_session_v3.py, mobile_v3.py, tsync_v3.py

CONFIGURATION (edit this block)
  IP_PREFIX = "10" : direct convention, phone ID 36 -> 10.153.247.36   (ID unchanged)
  IP_PREFIX = "192": legacy convention, phone ID 1  -> 192.168.164.101 (ID + 100)
  NETWORK is the 3rd IP octet; PORT is the Gobannos HTTP port.

EXAMPLES
  py connect_v3.py                       # prints the example address(es) and scans the subnet with nmap
  import connect_v3 as connect
  connect.get_adress(36)                 # '10.153.247.36'
  connect.get_adress(36, network=200)    # '10.153.200.36'
  connect.get_adresslist([36, 37])       # ['10.153.247.36', '10.153.247.37']
"""
import subprocess
from pprint import pprint

# ===================== USER CONFIGURATION =====================
# IP_PREFIX = "192"       # legacy convention
IP_PREFIX = "10"          # direct convention (current set-up)
LEGACY_SECOND_OCTET = 168
TEN_SECOND_OCTET = 153
NETWORK = 164 if IP_PREFIX == "192" else 247
PORT = 8080
DEFAULT_PHONELIST = [101, 102] if IP_PREFIX == "192" else [36]
# ==============================================================

if IP_PREFIX == "192":
    LOCAL_PREFIX = f"192.{LEGACY_SECOND_OCTET}"
    PHONE_OFFSET = 100
elif IP_PREFIX == "10":
    LOCAL_PREFIX = f"10.{TEN_SECOND_OCTET}"
    PHONE_OFFSET = 0
else:
    raise ValueError("IP_PREFIX must be '192' or '10'.")


def ipbase(network=NETWORK):
    """First three address octets plus a trailing dot."""
    return f"{LOCAL_PREFIX}.{network}."


def get_adress(phone, network=NETWORK):
    """Phonefleet ID -> IP, using the selected convention."""
    last_octet = int(phone) + PHONE_OFFSET
    if not 1 <= last_octet <= 254:
        raise ValueError(f"Phone ID {phone} maps to invalid final octet {last_octet}.")
    return f"{ipbase(network)}{last_octet}"


def get_adresslist(phonelist, network=NETWORK):
    """List of phone IDs -> list of IPs."""
    return [get_adress(phone, network) for phone in phonelist]


def basefolder():
    """Retained for compatibility with older Phonefleet code."""
    return "Bic25/"


def scan(network=NETWORK):
    """Scan the configured /24 subnet with nmap and return the Phonefleet IDs found."""
    subnet = f"{LOCAL_PREFIX}.{network}.0/24"
    try:
        result = subprocess.run(["nmap", "-sn", subnet], text=True, capture_output=True, check=False)
    except FileNotFoundError:
        print("nmap is not installed or not on PATH.")
        return []
    if result.returncode:
        print("nmap scan failed.")
        if result.stderr:
            print(result.stderr.strip())
        return []
    pprint(result.stdout.splitlines())
    prefix = f"{LOCAL_PREFIX}.{network}."
    phones = []
    for line in result.stdout.splitlines():
        marker = "Nmap scan report for "
        if marker not in line:
            continue
        address = line.split(marker, 1)[1].strip()
        if address.endswith(")") and "(" in address:       # "hostname (IP)"
            address = address.rsplit("(", 1)[1][:-1]
        if not address.startswith(prefix):
            continue
        try:
            phone_id = int(address[len(prefix):]) - PHONE_OFFSET
            if phone_id > 0:
                phones.append(phone_id)
        except ValueError:
            print(f"Could not parse nmap address: {address}")
    return sorted(set(phones))


def connect():
    """Scan and return (addresses, phone IDs)."""
    phones = scan()
    addresses = get_adresslist(phones)
    print(phones)
    print(f"Number of phones connected: {len(phones)}")
    return addresses, phones


def test_connect():
    """Print ADB devices, if the Android Debug Bridge is installed."""
    result = subprocess.run(["adb", "devices"], text=True, capture_output=True, check=False)
    print(result.stdout)
    if result.stderr:
        print(result.stderr)


if __name__ == "__main__":
    print(f"IP scheme: {IP_PREFIX}; network: {NETWORK}")
    print("Example address(es):", get_adresslist(DEFAULT_PHONELIST))
    print("Discovered phone IDs:", scan())