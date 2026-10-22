Architecture for controlling Chipiron



* A smartphone (now FP5) on Gobannos. 

&nbsp;	Accessible from a local network knowing its ip (broadcast wifi ?)

&nbsp;	Accessible from anywhere using ZeroTier/ Qrvpn (need MAC adress for the correspondance)

&nbsp;	Can be ssh, required Termux running + phone identifier

* USB hub to

&nbsp;	Connect to an Arduino Mega2560

&nbsp;	Power the FP externally (external battery ?)

* Arduino Mega2560 with

&nbsp;	- DO 13, PWM, duty cycle from 10% (0 speed) to 90% (max speed)

&nbsp;	- AI1 : speed data from the controller

&nbsp;	- AI2 : current data from the controller

* Controller ESCON 70/10

&nbsp;	- Power 24V (not very sensitive to the value : should be good without Tension control)

&nbsp;	- Wired to EC Motor

&nbsp;	- Hall sensor connected (Current measurement -> power ? Depends on the input power ? **check**)

&nbsp;	- Optical encoder for speed (motor speed without the gear !)

&nbsp;	- Digital IN/OUT : can be configured from the ESCON software (configured on Babasse)

 	- Analog IN/OUT : can be configured from the ESCON software (configured on Babasse)/ AO2 : enable (High), AO1 : PWM is now in speed control mode

 	- Controlling options, in particular speed range, and output ranges from the ESCON Software -> need some test in a granular medium

* Motor

&nbsp;	- Spécifications available in the datasheet. Max Current 9A

&nbsp;	- Gear ratio : 1184/18 ~1/66





