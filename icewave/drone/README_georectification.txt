# How to geo-rectify an image obtained from an UAV ? 


## Required data and parameters : 
	- an image obtained from an UAV (Unmanned Aerial Vehicle) using a calibrated camera
	- f, focal length of the camera (in pixels)
	- camera sensor dimensions in pixels : Lx, Ly
	- arrays of pixel coordinates of the image or of any object that appears on the image
	- H, drone altitude
	- alpha_0, camera picth angle to the horizontal, 
	- latitude and longitude of the UAV
	- azimuth of the camera, between 0° and 360°, defined clockwise from the North direction

The camera focal length can be ccomputed thanks to the Camera Calibrator package of Matlab. 
Flight parameters associated to the drone can be extracted from the .SRT and .csv files as described in the file
README_drone_files.txt

## 1st : Backward projection - Metric coordinates 
We use the pinhole camera model, where the camera optic system is simplified to a camera pupil entrance and a sensor plane. 
All light rays reaching the camera sensor plane are assumed to go through the pupil entrance. 
This backward projection enables to convert any pixel coordinates (xp,yp) associated to the camera sensor to metric coordinates (X,Y) 
associated to the captured ice plane. (see publication of data https://doi.org/10.5194/egusphere-2025-3304 for more detailed explanations and schematic)

This is achieved by a single function in the python module drone_projection.py : projection_real_space



## 2nd : Georeferencing - GPS position
Coordinates (X,Y) of the metric system associated to the UAV can then be converted into GPS coordinates. 
This is achieved thanks to the function from the module drone_projection.py : LatLong_coords_from_referencepoint
This function computes the latitude and longitude of any point using its distance and orientation from a point of knwon GPS coordinates. 

The georeferencing step is achieved through the following steps : 
	1 - compute GPS coordinates (Lat0,Long0) of the metric coordinate system center previously used. This point of 
	coordinates (X,Y) = (0,0) is associated to the center of the camera sensor (pixel coordinates (xp,yp)). This is achieved thanks 
	to the known GPS coordinate of the UAV, its azimuth, altitude and the pitch angle of the camera

	2 - compute distance and orientation of any point observed on the image to this reference point of coordinate (Lat0,Long0).
	Make use of the function LatLong_coords_from_referencepoint to compute the GPS coordinates 






