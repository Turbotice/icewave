#%%
import matplotlib.pyplot as plt
from Phidget22.Phidget import *
from Phidget22.Devices.VoltageRatioInput import *
import time
import os
#%%
#%matplotlib qt
#%%

# --- Configuration ---
name_output_file = r"C:\Users\Vasco Zanchi\Desktop\CASCADE26\tests_interfacage_test_bresilien\20261016\tests_bresilens\PhidgetBridge\test1.txt"
save = True

TIME_OUT = 5000  # 5s timeout
DATA_INTERVAL = 50  # 50ms sample frequency
remove_interval = 20  # Intervalle de suppression des points (en secondes)
dt = 0.01  # Temps entre les mises à jour (plus petit pour fluidité)

channels2plot = [0]  # Canaux à tracer

# Coefficients pour chaque canal
A = [1, 1, 1, 1]
B = [0, 0, 0, 0]

# --- Initialisation du fichier ---
if save:
    if os.path.exists(name_output_file):
        print('Fichier existant, ajout des données...')
        data_file = open(name_output_file, 'a')
        data_file.write('\nnew data:\n')
    else:
        print('Création d un nouveau fichier...')
        data_file = open(name_output_file, 'w')
        data_file.write('time vr0 vr1 vr2 vr3\n')
    data_file.flush()

# --- Initialisation des canaux Phidget ---
voltageRatioInputs = []
for channel in range(4):
    vr = VoltageRatioInput()
    vr.output_value = 0
    vr.setDeviceSerialNumber(680648)
    vr.setChannel(channel)

    # Callback dynamique
    def make_callback(ch):
        def callback(self, voltageRatio):
            self.output_value = (voltageRatio - B[ch]) / A[ch]
        return callback

    vr.setOnVoltageRatioChangeHandler(make_callback(channel))

    # ⚠️ OUVRIR LE CANAL AVANT DE LE CONFIGURER ⚠️
    vr.openWaitForAttachment(TIME_OUT)  # ✅ Attend que le périphérique soit attaché

    # ⚠️ CONFIGURATION MATÉRIELLE (APRÈS openWaitForAttachment) ⚠️
    vr.setBridgeGain(BridgeGain.BRIDGE_GAIN_128)  # ✅ OK maintenant
    vr.setDataInterval(DATA_INTERVAL)             # ✅ OK maintenant

    voltageRatioInputs.append(vr)

# --- Initialisation de la figure ---
#plt.ion()  # Mode interactif ON
fig, ax = plt.subplots()
ax.set_xlabel('Time [s]')
ax.set_ylabel('Voltage ratio')
ax.set_xlim(0, remove_interval)

# Initialisation des lignes (une par canal)
lines = []
for i, channel in enumerate(channels2plot):
    line, = ax.plot([], [], '.', color=f'tab:{["blue", "orange", "green", "purple"][i]}', label=f'vr{channel}')
    lines.append(line)
ax.legend()

# Stockage des données pour chaque ligne
times = [[] for _ in channels2plot]
values = [[] for _ in channels2plot]

# --- Boucle principale ---
try:
    t0 = time.time()
    while True:
        curr_time = time.time()
        t = curr_time - t0

        # Mise à jour des données
        for i, channel in enumerate(channels2plot):
            vr = voltageRatioInputs[channel]
            times[i].append(t)
            values[i].append(vr.output_value)

            # Mise à jour de la ligne correspondante
            lines[i].set_data(times[i], values[i])

            # Suppression des points trop anciens
            while len(times[i]) > 0 and t - times[i][0] > remove_interval:
                times[i].pop(0)
                values[i].pop(0)

        # Mise à jour des axes
        ax.set_xlim(max(0, t - remove_interval), t)
        ax.relim()  # Recalcule les limites des axes
        ax.autoscale_view()  # Ajuste automatiquement l'échelle

        # Écriture dans le fichier
        if save:
            voltages = [vr.output_value for vr in voltageRatioInputs]
            data_file.write(f"{t} {' '.join(map(str, voltages))}\n")
            data_file.flush()

        # Mise à jour de la figure
        fig.canvas.draw()
        fig.canvas.flush_events()
        time.sleep(dt)  # Petit délai pour éviter la surcharge CPU

finally:
    # Fermeture propre
    if save:
        data_file.close()
    for vr in voltageRatioInputs:
        vr.close()
    plt.ioff()  # Mode interactif OFF
    plt.show()  # Affiche la figure finale (optionnel)
# %%
