


def correspond_samplenum_acqnum(acq=1,serie=1):
    if serie!=None:
        fd = ((serie-1) * 3 + 1) + acq//3
    else:
        fd = (acq)//3
    sd = acq%3
    if sd==0:
        sd=3
    if fd==0:
        fd = 1
    return fd, sd