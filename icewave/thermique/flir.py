from flirimageextractor import FlirImageExtractor
import numpy as np
import glob
import matplotlib.pyplot as plt

flir = FlirImageExtractor()

filelist = glob.glob('exemples/*.jpg')
filename = filelist[0]

flir.process_image(flir_img_file=filename)
temperature = flir.get_thermal_np()  # Matrice des températures
max_temp = np.amax(temperature)
min_temp = np.amin(temperature)

plt.imshow(temperature, cmap='hot', interpolation='nearest')
plt.colorbar(label='Temperature (°C)')
plt.show()