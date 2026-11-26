"""
rw_data_v3.py -- tiny CSV read/write helpers (used to save the clock-synchronisation results).

DEPENDENCIES   stdlib only.
USED BY        tsync_v3.py, mobile_v3.py

FORMATS
  write_csv(file, rows)           rows = list of lists -> plain CSV (folder created if missing)
  writedict_csv(file, data)       data = {key: {col: value}} -> first column '#' = key, then one column per field
  read_csv(file) / csv2dict(table) read back (csv2dict converts numbers; use verbose=True to print the keys)

EXAMPLE
  import rw_data_v3 as rw
  rw.writedict_csv("tsync_1.txt", {36: {"tlag": 0.0123, "tstd": 1e-4}})
  table = rw.read_csv("tsync_1.txt");  data = rw.csv2dict(table)   # {36: {'tlag': 0.0123, 'tstd': 0.0001}}
"""
import csv
import os


def write_csv(filename, data):
    """Write a list of rows as CSV; creates the folder if needed."""
    folder = os.path.dirname(filename)
    if folder and not os.path.exists(folder):
        print(f"Creating folder {folder}")
        os.makedirs(folder)
    with open(filename, "w", newline="") as csvfile:
        writer = csv.writer(csvfile, delimiter=",", quotechar="|")
        for d in data:
            writer.writerow(d)


def writedict_csv(filename, data, symbol="#"):
    """Write {key: {col: value}} as CSV (header = symbol + the columns of the first entry). Empty dict: nothing."""
    if not data:
        print(f"writedict_csv: nothing to write in {filename}")
        return
    folder = os.path.dirname(filename)
    if folder:
        os.makedirs(folder, exist_ok=True)
    with open(filename, "w", newline="") as csvfile:
        writer = csv.writer(csvfile, delimiter=",", quotechar="|")
        first = next(iter(data.values()))
        writer.writerow([symbol] + list(first.keys()))
        for key, entry in data.items():
            writer.writerow([key] + [entry[k] for k in first.keys()])


def read_csv(filename, delimiter=","):
    """Read a CSV into a list of rows (strings)."""
    with open(filename, "r") as csvfile:
        return [row for row in csv.reader(csvfile, delimiter=delimiter, quotechar="|")]


def csv2dict(table, headerindex=0, symbol="#", verbose=False):
    """Rows -> dict. If the first cell is `symbol`: {key: {col: value}} with numbers converted; else {col: [values]}."""
    data = {}
    if table[0][0] == symbol:
        keys = table[0][1:]
        if verbose:
            print(keys)
        for tab in table[1:]:
            try:
                tab[0] = int(tab[0])
            except ValueError:
                pass
            data[tab[0]] = {}
            for t, key in zip(tab[1:], keys):
                if len(t.split(".")) > 2:
                    data[tab[0]][key] = t
                elif len(t.split(".")) == 2:
                    try:
                        data[tab[0]][key] = float(t)
                    except ValueError:
                        data[tab[0]][key] = t
                else:
                    try:
                        data[tab[0]][key] = int(t)
                    except ValueError:
                        data[tab[0]][key] = str(t)
    else:
        header = table[headerindex]
        data = {key: [] for key in header}
        for i in range(headerindex + 1, len(table)):
            for j, key in enumerate(header):
                data[key].append(table[i][j])
    return data