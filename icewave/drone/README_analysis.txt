# How to load and analyze displacement field computed by PIV algorithm ? 


## Results from PIV and post-processing 
	After performing PIV and post-processing, a file '*scaled.mat' should be available. It contains a matlab structure (equivalent to a python dictionnary), 
	with the following fields : 
	- Vx,Vy, (Vz for oblique views): displacement fields 
	- x,y,X,Y: spacial coordinates of the PIV boxes
	- t, UTC_t : time and UTC time
	- units : all units of the previous variables
	- PIXEL : pixel positions of PIV boxes on the camera sensor
	- SCALE : scale and factor scaling 
	- DRONE : drone parameters, height, angle, focal length, sensor size 
	- GPS : drone GPS position


## How to read it ?
	Matlab : load the '.mat' file as usual
	Python : use package 'matlab2python' to load the matlab structure and save it in a python structure. 
	An example is given in the python script : 'example_redressement_PIV_oblique.py'



## Useful functions for processing 
	Several functions are developped to analyze displacement fields both in matlab and python. 

	In Matlab, these functions are contained in the folder 'Drone_banquise_analysis', a few example are described below :
	- get_histogram.m : histogram of mean displacement, enable to check validity of PIV analysis
	- supress_quadratic_noise.m : supress noise due to drone motion. We fit each frame by a quadratic field and suppress
		this quadratic field. 
	- movie_velocity_field.m : generate a movie of the velocity field 
	- temporal_FFT.m : perform Fourier transform along time dimension
	- space_time_spectrum_Efk.m : perform Fourier transform in time and space, radial average in space dimension.
		Can be useful to plot dispersion relation

	and many others... 

	In Python, most of these functions are contained in two different python modules : 
	- drone_projection.py : contains several functions, mostly to perform conversion between different coordinates system, 
	of the camera sensor or the real space. 
	- Fourier_tools.py : contains several functions to perform Fourier transform and other things. In this module you can 
	find the following functions : 
		+ temporal_FFT, perform Fourier transform along time dimension
		+ supress_quadratic_noise, supress noise due to drone motion. We fit each frame by a quadratic field and suppress
		this quadratic field.
		+ space_time_spectrum, perform Fourier transform in time and space, radial average in space dimension.
		Can be useful to plot dispersion relation


