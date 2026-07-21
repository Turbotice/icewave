import numpy as np

def detrend_along_time_axis(v, window_size=120): # à peu près nombre de frames par periode d'onde):
    window = np.ones(window_size) / window_size
    v_detrended = np.zeros_like(v)

    # Supposons vz.shape = (N_x, N_y, N_t)
    for i in range(v.shape[0]):
        for j in range(v.shape[1]):
            v_detrended[i, j, :] = v[i, j, :] - np.convolve(v[i, j, :], window, mode='same')
    return v_detrended

vx_detrended = detrend_along_time_axis(vx)
vy_detrended = detrend_along_time_axis(vy)
vz_detrended = detrend_along_time_axis(vz)


ux = np.cumsum(vx_detrended, axis=2)*dt
uy = np.cumsum(vy_detrended, axis=2)*dt
uz = np.cumsum(vz_detrended, axis=2)*dt

%matplotlib qt
plt.figure(figsize=(15,10))
plt.plot(uz[10,10,:], '-b', alpha=0.5, label='numerical integration of detrended velocity field')
plt.plot(uz_trend[10,10,:], '-r', alpha=0.5, label='numerical integration of raw data')
plt.plot(uz[12,12,:], '-b', alpha=0.5)
plt.plot(uz_trend[12,12,:], '-r', alpha=0.5)
plt.legend()
plt.show()



