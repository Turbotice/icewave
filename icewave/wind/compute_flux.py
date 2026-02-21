import numpy as np
import xarray as xr


def compute_mean_X(file_nc, 
                   X='U', 
                   period=600,
                   sampling=5,
                   method='block',
                   update_nc=True):
    """
    Computes the average of the variables in ds, with a period of 10 min by defaults.
    period in secondes.
    Returns the dataset with a variable 'mean_X' added.
    """
    
    ds = xr.open_dataset(file_nc)
    length_window = 10 * int(period / sampling)
    # Check if work has been done
    if 'mean_'+X in ds.keys():
        print(f'Mean for {X} already in the file, skipping ...')
        return
    
    
    
    if method=='block':
        ds['mean_'+X] = ds[X]*0.
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

    elif method=='moving':
        length_window = int(period * sampling)
        ds['mean_'+X] = ds[X].rolling(time=length_window, center=True).mean('time')


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


def rotation(file_nc, mean_kwargs={}, update_nc=True):
    """
    Rotation of the reference frame

    INPUTS:
        file_nc: str, the netcdf file to update
        update_nc: bool, force overwrite

    Goal is to have:
        <V> = <W> = <v'w'> = 0
    So we need:
    - a rotation around z1 axis (eta) => <V> = 0
    - a rotation around y1 axis (theta) => <W> = 0
    

    We use the convention by Lee et al. [1]:
        Intrument reference frame: U1
        Intermediate frame <V>=<W>=0: U2
        Natural wind frame <V>=<W>=<w'v'>=0: U
    
    REFERENCES:

    [1] Lee, X., Massman, W., & Law, B. (Eds.). (2005). 
        Handbook of Micrometeorology (Vol. 29). Springer Netherlands. 
        https://doi.org/10.1007/1-4020-2265-4

    """
    ds = xr.open_dataset(file_nc)

    # precompute some values
    CE = (ds.mean_U1 / np.sqrt(ds.mean_U1**2+ds.mean_V1**2)).values
    SE = (ds.mean_V1 / np.sqrt(ds.mean_U1**2+ds.mean_V1**2)).values
    CT = (np.sqrt(ds.mean_U1**2+ds.mean_V1**2)
          /
          np.sqrt(ds.mean_U1**2+ds.mean_V1**2+ds.mean_W1**2)).values
    ST = ( ds.mean_W1
            /
          np.sqrt(ds.mean_U1**2+ds.mean_V1**2+ds.mean_W1**2)).values
        
    # Rotation <V>=<W>=0 in one go
    ds['U2'] = ds.U1*CT*CE +
                ds.V1*CT*SE +
                ds.W1*ST
    ds['V2'] = ds.V1*CE - ds.U1*SE
    ds['W2'] = ds.W1*CT - ds.U1*ST*CE - ds.V1*ST*SE

    # Rotation for <u'w'>=0
    compute_mean_X(file_nc, 
                X=var, 
                period=avg_period, 
                sampling=sampling,
                method=avg_method,
                update_nc=True)
    mean_w2v2 = 0
    beta = 0.5*np.tan()

    raise Exception('WIPPPP')

    # Saving 
    if update_nc:
        ds.load()
        ds.close() # <- required, else cant append to file
        ds.to_netcdf(file_nc, mode="a")
















