import numpy as np
import xarray as xr


def compute_mean_X(file_nc, 
                   X='U', 
                   period=600, 
                   method='block',
                   update_nc=True):
    """
    Computes the average of the variables in ds, with a period of 10 min by defaults.
    period in secondes.
    Returns the dataset with a variable 'mean_X' added.
    """
    
    ds = xr.open_dataset(file_nc)

    # Check if work has been done
    if 'mean_'+X in ds.keys():
        print(f'Mean for {X} already in the file, skipping ...')
        return
    
    ds['mean_'+X] = ds[X]*0.
    time = ds.coords['time']
    
    if method=='block':
        t0 = ds.time[0]
        tlast = ds.time[-1]
        delta =  np.timedelta64(period, 's')
        t1 = t0 + delta
        while t0 < tlast:
            # select only 1 block
            dsm = ds.sel(time=slice(t0,t1))
            # compute mean
            meanX = dsm[X].mean('time').values
            # assign value
            ds["mean_"+X].loc[dict(time=slice(t0, t1))] = meanX

            # next iteration
            t0 = t1
            t1 = t1 + delta  
    else:
        raise Exception(f'This method {method} isnt coded') 


    # Saving
    if update_nc:
        ds.load()
        ds.close() # <- required, else cant append to file
        ds.to_netcdf(file_nc, mode="a")


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
    liste_variables = ['U','V','W','T']
    for var in liste_variables:
        ds1 = function(ds, X=var, period=period)

    return ds1



