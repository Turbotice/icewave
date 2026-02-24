import numpy as np
import xarray as xr


def mean_operator(ds, X='U', 
                   period=600,
                   sampling=5,
                  method='block'):
    """
    Computes the average of the variables in ds, with a period of 10 min by defaults.
    period in secondes.

    INPUTS:
        ds: xr.Dataset
        X: str, variable to average
        sampling: float, sampling rate of the instrument
        method: str, what kind of average to use
    OUPUTS:
        xr.DataArray of the average variable
    """
    # Check if work has been done
    if 'mean_'+X in ds.keys():
        print(f'Mean for {X} already in the file, skipping ...')
        return ds['mean_'+X]

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
        raise Exception(f'This method ({method}) isnt coded')
    return ds['mean_'+X]





def rotation_operator(ds, mean_kwargs={}):
    """
    Rotation of the reference frame

    INPUTS:
        ds: xr.Dataset 
        mean_kwargs: dict, passed to mean_operator
    OUTPUTS:
        ds: xr.Dataset, updated dataset with rotation applied
    Goal is to have:
        <V> = <W> = <v'w'> = 0
    So we need:
    - a rotation around z1 axis (eta) => <V> = 0
    - a rotation around y1 axis (theta) => <W> = 0
    

    We use the convention by Lee et al. [1]:
        Intrument reference frame: U1
        Intermediate frame <V>=<W>=0: U2
        Natural wind frame <V>=<W>=<w'v'>=0: U
    
    and: 
        <U>  <=>  mean_U
        u'   <=>  up

    REFERENCES:

    [1] Lee, X., Massman, W., & Law, B. (Eds.). (2005). 
        Handbook of Micrometeorology (Vol. 29). Springer Netherlands. 
        https://doi.org/10.1007/1-4020-2265-4

    """
    # first compute averages in instrument reference frame
    for var in ['U1','V1','W1']:
        if 'mean_'+var not in ds.keys():
            print(var+f' not in the dataset, computing mean_{var} ...')

            # ds['mean_'+var]
            ds1 = mean_operator(ds, 
                                            X=var,
                                            period=mean_kwargs['period'], 
                                            sampling=mean_kwargs['sampling'], 
                                            method=mean_kwargs['method'] )
    
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
    ds['U2'] = (ds.U1*CT*CE +
                ds.V1*CT*SE +
                ds.W1*ST )
    ds['V2'] = ds.V1*CE - ds.U1*SE
    ds['W2'] = ds.W1*CT - ds.U1*ST*CE - ds.V1*ST*SE

    # Rotation for <u'w'>=0
    for var in ['V2','W2']:
        ds['mean_'+var] = mean_operator(ds, 
                                        X=var,
                                        period=mean_kwargs['period'], 
                                        sampling=mean_kwargs['sampling'], 
                                        method=mean_kwargs['method'] )
    ds['vp2wp2'] = (ds.V2 - ds.mean_V2)*(ds.W2 - ds.mean_W2)
    ds['vp2vp2'] = (ds.V2 - ds.mean_V2)*(ds.V2 - ds.mean_V2)
    ds['wp2wp2'] = (ds.W2 - ds.mean_W2)*(ds.W2 - ds.mean_W2)

    for var in ['vp2wp2','vp2vp2','wp2wp2']:
         ds['mean_'+var] = mean_operator(ds, 
                                        X=var,
                                        period=mean_kwargs['period'], 
                                        sampling=mean_kwargs['sampling'], 
                                        method=mean_kwargs['method'] )
    # computing beta
    beta = 2*np.tan(2 * ds.mean_vp2wp2/ (ds.mean_vp2vp2 - ds.mean_wp2wp2))
    beta = 1/beta
    
    CB = np.cos(beta)
    SB = np.sin(beta)

    # Final values
    ds['U'] = ds['U2']
    ds['V'] = ds['V2']*CB + ds['W2']*SB
    ds['W'] = ds['W2']*CB - ds['V2']*SB

    # No rotation for T (scalar)
    if 'T' not in ds.keys():
        ds = ds.rename({'T1':'T'})

    return ds


def compute_flux_wx(ds, X='U', period=600):
    """
    Compute the flux w'x' with x in [U,V,W,T]
    Returns the dataset with a variable 'flux_wx' added.
    """
    raise Exception('to do')


    ds1 = ds
    if 'mean_'+X not in ds.keys():
        ds1 = compute_mean(ds, period)
    
    raise Exception('to do')

    return ds1




# ==================================
# TOOLBOX
# ==================================
#
def compute_and_save(file_nc, operator, update_nc=True, **kwargs):
    """
    Generic wrapper for any operator (mean_operator, rotation_operator, etc.)
    """
    ds = xr.open_dataset(file_nc)
    
    # Compute
    ds = operator(ds, **kwargs)
    
    # Save
    if update_nc:
        ds.load()
        ds.close()
        ds.to_netcdf(file_nc, mode="a")


       















