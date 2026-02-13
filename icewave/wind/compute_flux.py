import numpy as np
import xarray as xr


def compute_mean_X(ds, X='U', period=600):
    """
    Computes the average of the variables in ds, with a period of 10 min by defaults.
    period in secondes.
    Returns the dataset with a variable 'mean_X' added.
    """
    raise Exception('to do')


def compute_flux_wx(ds, X='U', period=600):
    """
    Compute the flux w'x' with x in [U,V,W,T]
    Returns the dataset with a variable 'flux_wx' added.
    """
    ds1 = ds
    if 'mean_'+X not in ds.keys():
        ds1 = compute_mean(ds, period)

    raise Exception('to do')

    return ds1

def compute_all_operator(ds, function, period=600):
    """
    Computes all f(X) quantities and returns an updated dataset
    """
    ds1 = ds
    liste_variables = ['S','U','V','W','T']
    for var in liste_variables:
        ds1 = function(ds, X=var, period)

    return ds1



