#%%
import matplotlib.pyplot as plt
from Phidget22.Phidget import *
from Phidget22.Devices.VoltageRatioInput import *
import time
import os
#%%
%matplotlib qt
#%%
####################################
name_output_file = "C:/Users/Vasco Zanchi/Desktop/CASCADE26/calibration_capteur_force/loadcell_n1080/discharge_loadcell_n1080_mass788g_3.txt"

save = False

def create_new_data_file(name_output_file):
    file = open(name_output_file,'w')
    file.write('time vr0 vr1 vr2 vr3\n')
    file.close()

def add_text_to_data_file(name_output_file, text):
    file = open(name_output_file,'a')
    file.write(text)
    file.close()

if (os.path.exists(name_output_file)) & (save==True):
    print('output file already exists, adding in that file..')
    file = open(name_output_file,'a')
    file.write('new data :\n')
    file.close()
elif (not(os.path.exists(name_output_file))) & (save==True):    
    create_new_data_file(name_output_file)
#####################################




TIME_OUT = 5000  # 5s before it throws a timeout exception
DATA_INTERVAL = 50  # 1000ms sample frequency

A0 = 1
B0 = 0

A1 = 1
B1 = 0

A2 = 1
B2 = 0

A3 = 1
B3 = 0

dt = 0.1 # number of seconds before updating the figure
remove_interval = 20 # in seconds

channels2plot = [0] # on schoisit les channels qui nous intéresse (par exemple si il y a un appareil seulement sur ch0 on met [0])

def onVoltageRatioChange0(self, voltageRatio): 
    OutputValue = (voltageRatio - (B0)) / (A0)
    self.output_value = OutputValue

def onVoltageRatioChange1(self, voltageRatio):
    OutputValue = (voltageRatio - (B1)) / (A1)
    self.output_value = OutputValue

def onVoltageRatioChange2(self, voltageRatio):
    OutputValue = (voltageRatio - (B2)) / (A2)
    self.output_value = OutputValue

def onVoltageRatioChange3(self, voltageRatio):
    OutputValue = (voltageRatio - (B3)) / (A3)
    self.output_value = OutputValue

def main():
    # Initialize VoltageRatioInput objects for the first bridge
    voltageRatioInput0 = VoltageRatioInput()
    voltageRatioInput0.output_value = 0
    voltageRatioInput0.setDeviceSerialNumber(680648)
    voltageRatioInput0.setChannel(0)
    voltageRatioInput0.setOnVoltageRatioChangeHandler(onVoltageRatioChange0)
    voltageRatioInput0.openWaitForAttachment(TIME_OUT)
    voltageRatioInput0.setBridgeGain(BridgeGain.BRIDGE_GAIN_128)
    voltageRatioInput0.setDataInterval(DATA_INTERVAL)

    voltageRatioInput1 = VoltageRatioInput()
    voltageRatioInput1.output_value = 0
    voltageRatioInput1.setDeviceSerialNumber(680648)
    voltageRatioInput1.setChannel(1)
    voltageRatioInput1.setOnVoltageRatioChangeHandler(onVoltageRatioChange1)
    voltageRatioInput1.openWaitForAttachment(TIME_OUT)
    voltageRatioInput1.setBridgeGain(BridgeGain.BRIDGE_GAIN_128)
    voltageRatioInput1.setDataInterval(DATA_INTERVAL)

    voltageRatioInput2 = VoltageRatioInput()
    voltageRatioInput2.output_value = 0
    voltageRatioInput2.setDeviceSerialNumber(680648)
    voltageRatioInput2.setChannel(2)
    voltageRatioInput2.setOnVoltageRatioChangeHandler(onVoltageRatioChange2)
    voltageRatioInput2.openWaitForAttachment(TIME_OUT)
    voltageRatioInput2.setBridgeGain(BridgeGain.BRIDGE_GAIN_128)
    voltageRatioInput2.setDataInterval(DATA_INTERVAL)

    voltageRatioInput3 = VoltageRatioInput()
    voltageRatioInput3.output_value = 0
    voltageRatioInput3.setDeviceSerialNumber(680648)
    voltageRatioInput3.setChannel(3)
    voltageRatioInput3.setOnVoltageRatioChangeHandler(onVoltageRatioChange3)
    voltageRatioInput3.openWaitForAttachment(TIME_OUT)
    voltageRatioInput3.setBridgeGain(BridgeGain.BRIDGE_GAIN_128)
    voltageRatioInput3.setDataInterval(DATA_INTERVAL)

    voltageRatioInputs_bridge1 = [voltageRatioInput0,voltageRatioInput1,voltageRatioInput2,voltageRatioInput3]

    try:
        t0 = time.time()
        plt.figure()
        plt.xlabel('Time [s]')
        plt.ylabel('Voltage ratio')
        while True:
            curr_time = time.time()
            t = curr_time - t0
            
            # Read and plot data from the first bridge
            V0 = voltageRatioInput0.output_value
            V1 = voltageRatioInput1.output_value
            V2 = voltageRatioInput2.output_value
            V3 = voltageRatioInput3.output_value
            if 0 in channels2plot:
                plt.plot(t, voltageRatioInput0.output_value, '.', color='tab:blue', label='voltageRatioInput0')
            if 1 in channels2plot:
                plt.plot(t, voltageRatioInput1.output_value, '.', color='tab:orange', label='voltageRatioInput1')
            if 2 in channels2plot:
                plt.plot(t, voltageRatioInput2.output_value, '.', color='tab:green', label='voltageRatioInput2')
            if 3 in channels2plot:
                plt.plot(t, voltageRatioInput3.output_value, '.', color='tab:purple', label='voltageRatioInput3')
           
            plt.xlim(t-remove_interval, t)
            plt.pause(dt)
            t += dt            
            plt.draw()

            # Remove old data points
            lines_to_remove = []
            lines = plt.gca().get_lines()
            print("number of points on the figure : ", len(lines))
            for line in lines:
                x_data = line.get_xdata()[0]
                if t - x_data > remove_interval:
                    lines_to_remove.append(line)
            for line in lines_to_remove:
                line.remove()
                
            if save: 
                # Save data to file
                voltages = [voltageRatioInput.output_value for voltageRatioInput in voltageRatioInputs_bridge1]
                add_text_to_data_file(name_output_file, f"{t} {' '.join(map(str, voltages))}\n")

    finally:
        # Close connections
        for voltageRatioInput in voltageRatioInputs_bridge1:
            voltageRatioInput.close()
        plt.show()

if __name__ == "__main__":
    main()
# %%
