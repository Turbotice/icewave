**How to use Gobannos to control the motor**



⦁	A smartphone (now FP4) on Gobannos. 

&nbsp;	Accessible from a local network knowing its ip (broadcast wifi ?)

&nbsp;	Accessible from anywhere using ZeroTier (need MAC adress for the correspondance)

&nbsp;	Can be ssh, required Termux running + phone identifier



Beware of the order you connect the elements :

* Launch Gobannos
* Wired USB hub
* Wired Arduino
* Add power supply



**Run acquisition :**

use mobile.py (in phonefleet) or in a browser

http://192.168.network.phoneip/status (check status)

http://192.168.network.phoneip/start (start acquisition)



Control the motor speed:

you can set the value of the PWM using an url request, the syntax is the following :

http://192.168.network.phoneip/usb-cmd/cXX



where XX is the value of the PWM, from 10 (0 speed) to 90 (max speed).
Max speed is currently set to 2640rpm before reduction, corresponding to ~40rpm

**Check values**



The output values of the speed and the current are given by the controller to the Arduino board, and then Streamed to Gobannos

Speed comes first, then Current after a delay of 50ms. The waiting time for the next loop is 150ms+

The time stamp should then look like

t0	speed

t0+50ms	current

t1=t1+~200ms	speed



If several values follows briefly, (<50ms) they correspond to the bits of the same output (Arduino DI are in 12 bits)

Would require postprocessing to reconstruct the value





