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
name_output_file = 'C:/Users/Vasco/test1thermocouples.txt'

save = False

def create_new_data_file(name_output_file):
    file = open(name_output_file,'w')
    file.write('time V0 V1 V2 V3 V4 V5 V6 V7\n')
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
else:    
    create_new_data_file(name_output_file)
#####################################




TIME_OUT = 5000  # 5s before it throws a timeout exception
DATA_INTERVAL = 50  # 1000ms sample frequency

A0 = 1
B0 = 0#-0.7e-5#-2.7e-5

A1 = 1
B1 = 0#-0.1e-5#-2.7e-5

A2 = 1
B2 = 0#-0.5e-5 #-3.7e-5

A3 = 1
B3 = 0#-3.4e-5

dt = 0.5 # number of seconds before updating the figure
remove_interval = 20 # in seconds

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


# Define functions for the second Phidget bridge channels
def onVoltageRatioChange4(self, voltageRatio): 
    OutputValue = (voltageRatio - (B0)) / (A0)
    self.output_value = OutputValue

def onVoltageRatioChange5(self, voltageRatio):
    OutputValue = (voltageRatio - (B1)) / (A1)
    self.output_value = OutputValue

def onVoltageRatioChange6(self, voltageRatio): 
    OutputValue = (voltageRatio - (B2)) / (A2)
    self.output_value = OutputValue

def onVoltageRatioChange7(self, voltageRatio):
    OutputValue = (voltageRatio - (B3)) / (A3)
    self.output_value = OutputValue

def main():
    # Initialize VoltageRatioInput objects for the first bridge
    voltageRatioInput0 = VoltageRatioInput()
    voltageRatioInput0.output_value = 0
    voltageRatioInput0.setDeviceSerialNumber(586043)
    voltageRatioInput0.setChannel(0)
    voltageRatioInput0.setOnVoltageRatioChangeHandler(onVoltageRatioChange0)
    voltageRatioInput0.openWaitForAttachment(TIME_OUT)
    voltageRatioInput0.setBridgeGain(BridgeGain.BRIDGE_GAIN_128)
    voltageRatioInput0.setDataInterval(DATA_INTERVAL)

    voltageRatioInput1 = VoltageRatioInput()
    voltageRatioInput1.output_value = 0
    voltageRatioInput1.setDeviceSerialNumber(586043)
    voltageRatioInput1.setChannel(1)
    voltageRatioInput1.setOnVoltageRatioChangeHandler(onVoltageRatioChange1)
    voltageRatioInput1.openWaitForAttachment(TIME_OUT)
    voltageRatioInput1.setBridgeGain(BridgeGain.BRIDGE_GAIN_128)
    voltageRatioInput1.setDataInterval(DATA_INTERVAL)

    voltageRatioInput2 = VoltageRatioInput()
    voltageRatioInput2.output_value = 0
    voltageRatioInput2.setDeviceSerialNumber(586043)
    voltageRatioInput2.setChannel(2)
    voltageRatioInput2.setOnVoltageRatioChangeHandler(onVoltageRatioChange2)
    voltageRatioInput2.openWaitForAttachment(TIME_OUT)
    voltageRatioInput2.setBridgeGain(BridgeGain.BRIDGE_GAIN_128)
    voltageRatioInput2.setDataInterval(DATA_INTERVAL)

    voltageRatioInput3 = VoltageRatioInput()
    voltageRatioInput3.output_value = 0
    voltageRatioInput3.setDeviceSerialNumber(586043)
    voltageRatioInput3.setChannel(3)
    voltageRatioInput3.setOnVoltageRatioChangeHandler(onVoltageRatioChange3)
    voltageRatioInput3.openWaitForAttachment(TIME_OUT)
    voltageRatioInput3.setBridgeGain(BridgeGain.BRIDGE_GAIN_128)
    voltageRatioInput3.setDataInterval(DATA_INTERVAL)

    voltageRatioInput4 = VoltageRatioInput()
    voltageRatioInput4.output_value = 0
    voltageRatioInput4.setDeviceSerialNumber(680648)
    voltageRatioInput4.setChannel(0)
    voltageRatioInput4.setOnVoltageRatioChangeHandler(onVoltageRatioChange0)
    voltageRatioInput4.openWaitForAttachment(TIME_OUT)
    voltageRatioInput4.setBridgeGain(BridgeGain.BRIDGE_GAIN_128)
    voltageRatioInput4.setDataInterval(DATA_INTERVAL)

    voltageRatioInput5 = VoltageRatioInput()
    voltageRatioInput5.output_value = 0
    voltageRatioInput5.setDeviceSerialNumber(680648)
    voltageRatioInput5.setChannel(1)
    voltageRatioInput5.setOnVoltageRatioChangeHandler(onVoltageRatioChange1)
    voltageRatioInput5.openWaitForAttachment(TIME_OUT)
    voltageRatioInput5.setBridgeGain(BridgeGain.BRIDGE_GAIN_128)
    voltageRatioInput5.setDataInterval(DATA_INTERVAL)

    voltageRatioInput6 = VoltageRatioInput()
    voltageRatioInput6.output_value = 0
    voltageRatioInput6.setDeviceSerialNumber(680648)
    voltageRatioInput6.setChannel(2)
    voltageRatioInput6.setOnVoltageRatioChangeHandler(onVoltageRatioChange2)
    voltageRatioInput6.openWaitForAttachment(TIME_OUT)
    voltageRatioInput6.setBridgeGain(BridgeGain.BRIDGE_GAIN_128)
    voltageRatioInput6.setDataInterval(DATA_INTERVAL)

    voltageRatioInput7 = VoltageRatioInput()
    voltageRatioInput7.output_value = 0
    voltageRatioInput7.setDeviceSerialNumber(680648)
    voltageRatioInput7.setChannel(3)
    voltageRatioInput7.setOnVoltageRatioChangeHandler(onVoltageRatioChange3)
    voltageRatioInput7.openWaitForAttachment(TIME_OUT)
    voltageRatioInput7.setBridgeGain(BridgeGain.BRIDGE_GAIN_128)
    voltageRatioInput7.setDataInterval(DATA_INTERVAL)

    voltageRatioInputs_bridge1 = [voltageRatioInput0,voltageRatioInput1,voltageRatioInput2,voltageRatioInput3]
    voltageRatioInputs_bridge2 = [voltageRatioInput4,voltageRatioInput5,voltageRatioInput6,voltageRatioInput7]

    try:
        t0 = time.time()
        plt.figure()
        plt.xlabel('Time (s)')
        plt.ylabel('Voltage')
        while True:
            curr_time = time.time()
            t = curr_time - t0
            
            # Read and plot data from the first bridge
            V0 = voltageRatioInput0.output_value
            V1 = voltageRatioInput1.output_value
            V2 = voltageRatioInput2.output_value
            V3 = voltageRatioInput3.output_value
            V4 = voltageRatioInput4.output_value
            V5 = voltageRatioInput5.output_value
            V6 = voltageRatioInput6.output_value
            V7 = voltageRatioInput7.output_value
            
            plt.plot(t, voltageRatioInput0.output_value, '.', color='tab:blue', label='voltageRatioInput0')
            plt.plot(t, voltageRatioInput1.output_value, '.', color='tab:orange', label='voltageRatioInput1')
            plt.plot(t, voltageRatioInput2.output_value, '.', color='tab:green', label='voltageRatioInput2')
            plt.plot(t, voltageRatioInput3.output_value, '.', color='tab:purple', label='voltageRatioInput3')
            plt.plot(t, voltageRatioInput4.output_value, '^', color='tab:blue', label='voltageRatioInput0')
            plt.plot(t, voltageRatioInput5.output_value, '^', color='tab:orange', label='voltageRatioInput1')
            plt.plot(t, voltageRatioInput6.output_value, '^', color='tab:green', label='voltageRatioInput2')
            plt.plot(t, voltageRatioInput7.output_value, '^', color='tab:purple', label='voltageRatioInput3')
                
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
                voltages = [voltageRatioInput.output_value for voltageRatioInput in voltageRatioInputs_bridge1] + \
                           [voltageRatioInput.output_value for voltageRatioInput in voltageRatioInputs_bridge2]
                add_text_to_data_file(name_output_file, f"{t} {' '.join(map(str, voltages))}\n")

    finally:
        # Close connections
        for voltageRatioInput in voltageRatioInputs_bridge1:
            voltageRatioInput.close()
        for voltageRatioInput in voltageRatioInputs_bridge2:
            voltageRatioInput.close()
        plt.show()

if __name__ == "__main__":
    main()

# %%
