#%%
from djilogparser import DJIParser

# Path to your DJI TXT flight record
log_file = r"C:\Users\Vasco Zanchi\Desktop\DJIFlightRecord_2026-02-14_[18-39-10].txt"

parser = DJIParser(log_file)
data = parser.parse()

# Print available fields
print(data.keys())

# Example: print GPS records
for record in data["gps"]:
    print(record)
#%%
from djilogparser import DJIParser
import csv

log_file = r"C:\Users\Vasco Zanchi\Desktop\DJIFlightRecord_2026-02-14_[18-39-10].txt"
parser = DJIParser(log_file)
data = parser.parse()

with open("gps_output.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["latitude", "longitude", "altitude", "time"])

    for r in data["gps"]:
        writer.writerow([r["lat"], r["lon"], r["alt"], r["timestamp"]])

print("Export complete")
