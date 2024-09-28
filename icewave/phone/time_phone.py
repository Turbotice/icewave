import datetime


def get_time(t0):
    t00, t01 = t0.split(" UTC")[0].split(".")
    date = datetime.datetime.strptime(t00.split(" UTC")[0], "%Y-%m-%d %H:%M:%S")
    tphone = date.timestamp() + int(t01)/1000-3600
    ts = datetime.datetime.fromtimestamp(tphone)
    return tphone,ts
