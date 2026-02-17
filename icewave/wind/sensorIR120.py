import numpy as np
import matplotlib.pyplot as plt

def T_t(R): 
    """
    Thermistor temperature (Kelvin)
    """

    # Calibration constants
    A = 9.325752E-4
    B = 2.212819E-4
    C = 1.333627E-7
    
    return 1/(A+B*np.log(R)+C*np.log(R)**3)


R = np.arange(50, 200)
fig, ax = plt.subplots(1,1,figsize = (3,3), constrained_layout=True, dpi=100)
ax.plot(R, T_t(R*1000))
ax.hlines(273, 0, 1000, color='gray', ls='--')
ax.set_ylabel('T (K)')
ax.set_xlabel('R (kOhms)')
ax.set_xlim([50,200])
plt.show()
