"""
This file defines the different sensor types.

"""


class Anemometer:
    def __init__(self, fs, files, name, description, notes):
        self.fs = fs        # sampling rate
        self.files = files
        self.name = name
        self.description = description
        self.notes = notes

class T_sensor:
    def __init__(self, fs, files, name, kind):
        self.fs = fs
        self.fs = files
        self.name = name
        self.kind = kind
