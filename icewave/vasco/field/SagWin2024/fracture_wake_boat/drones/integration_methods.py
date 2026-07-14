#%%

window_size=120
window = np.ones(window_size) / window_size
vz_detrended = np.zeros_like(vz)

# Supposons vz.shape = (N_x, N_y, N_t)
for i in range(vz.shape[0]):
    for j in range(vz.shape[1]):
        vz_detrended[i, j, :] = vz[i, j, :] - np.convolve(vz[i, j, :], window, mode='same')

dt = 1/dict_stereo_pivdata['SCALE']['facq_t']

uz = np.cumsum(vz_detrended, axis=2) * dt

uz_trend = np.cumsum(vz, axis=2) * dt

%matplotlib qt
plt.figure(figsize=(15,10))
plt.plot(uz[10,10,:], '-b', alpha=0.5, label='numerical integration of detrended velocity field')
plt.plot(uz_trend[10,10,:], '-r', alpha=0.5, label='numerical integration of raw data')
plt.plot(uz[12,12,:], '-b', alpha=0.5)
plt.plot(uz_trend[12,12,:], '-r', alpha=0.5)
plt.legend()
plt.show()



