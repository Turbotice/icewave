# %% import libraries
import numpy as np
import matplotlib.pyplot as plt

# %% Load data
a = np.loadtxt('C:/Users/Vasco/OneDrive - Université de Paris/Bureau/stage_M2/profil_a.txt')
b = np.loadtxt('C:/Users/Vasco/OneDrive - Université de Paris/Bureau/stage_M2/profil_b.txt')

# %% Plot data
plt.plot(a)
b = np.hstack((a[10:],a[:10]))
plt.plot(b)
plt.show()

# %% Parameters
Vmax = 16

# %% Initialize arrays
P1 = np.zeros(2*Vmax-1)
P2 = np.zeros(2*Vmax-1)
print(len(P1))
# %% Calculate cross-correlation
for v in range(-Vmax+1, Vmax):
    P1[v+Vmax-1] = np.sum((a[Vmax-v-1:len(a)-Vmax-v+1] - b[Vmax-1:len(a)-Vmax+1])**2)
    P2[v+Vmax-1] = np.sum((b[Vmax-v-1:len(a)-Vmax-v+1] - a[Vmax-1:len(a)-Vmax+1])**2)

plt.plot(a[Vmax-v:len(a)-Vmax-v+1])
plt.plot(b[Vmax-v:len(a)-Vmax-v+1])
P2 = np.flip(P2)
P = P1 + P2

# %% Plot cross-correlation
plt.plot(np.arange(1, len(P)+1)-Vmax, P, '+')
plt.show()

# %% Detection sub-pixel
i_min = np.argmin(P)

delta_i = (P[i_min+1] - P[i_min-1]) / (2 * (2 * P[i_min] - P[i_min+1] - P[i_min-1]))

# %% Lag value obtained
lag = i_min - Vmax + delta_i + 1 # il faut rajouter 1 car argmin donne l'indice du minimum mais en convention python

print("Lag:", lag)


# %% define a function that detects the wave front with subpixelar precision
def detection_front_sub_pix(p_a,p_b,Vmax=16):
    P1 = np.zeros(2*Vmax-1)
    P2 = np.zeros(2*Vmax-1)
    for v in range(-Vmax+1, Vmax):
        P1[v+Vmax-1] = np.sum((a[Vmax-v-1:len(a)-Vmax-v+1] - b[Vmax-1:len(a)-Vmax+1])**2)
        P2[v+Vmax-1] = np.sum((b[Vmax-v-1:len(a)-Vmax-v+1] - a[Vmax-1:len(a)-Vmax+1])**2)
    P2 = np.flip(P2)
    P = P1 + P2

    if plot==True:
        plt.figure()
        plt.plot(a[Vmax-v:len(a)-Vmax-v+1])
        plt.plot(b[Vmax-v:len(a)-Vmax-v+1])
        plt.show()
        plt.plot(np.arange(1, len(P)+1)-Vmax, P, '+')
        plt.show()

    i_min = np.argmin(P)
    delta_i = (P[i_min+1] - P[i_min-1]) / (2 * (2 * P[i_min] - P[i_min+1] - P[i_min-1]))
    lag = i_min - Vmax + delta_i + 1 # il faut rajouter 1 car argmin donne l'indice du minimum mais en convention python

    return lag