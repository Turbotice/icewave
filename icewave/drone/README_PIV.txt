# How to use codes Drones/PIV_processing ?


This file describes different codes that are needed to perform a Digital Image Correlation (DIC) method on a video recorded with an UAV. 
The DIC algorithm is performed through the PIVlab software of Matlab. 


## From a video '.MP4' to a folder of images '.tiff' :
--> convert_mutlivideo2tiff.py : 
	converts a serie of video into a folder of images. Images are named by index, starting from 0. 

## Processing images using PIVlab : 

  ### Before starting : 
  --> DIC method needs several parameters : 
	Dt: time step between two images that are compared
	W: window size, in pixels, within which pixels of image i are compared to the corresponding pixels of image i+Dt
	b: time step between two image comparison : u[i] = img[i+1] - img[i] , u[i+1] = img[i+b+1] - img[i+b]


  --> Rules to follow : 
  The computed displacement field must respect the following conditions : 
	- u > 0.1 pix/frame => sufficiently high signal/noise ratio 
	- u < W/4 => pixel displacements must not be too high compared to the window size


  ### How to determine time step between two images ? 
  --> Open PIVlabe, import a serie of images (maximum 20 images) which are successively separated by [1,2,3,4,5,6,7...] images. Execute the PIVlab analysis with
	the desired parameters, pre-processing parameters, etc...
  --> Have a look at the obtained velocity fields for the different time steps. 
  --> Choose the value of Dt such that the computed fields check the rules previously mentionned
  --> This step can be repeated several times, with different sequences of images in order to define the best time step. 



 ### Perform PIV
 --> PIV_processing/automaticcycle_banquise_drone.m
	Enables to perform DIC algorithm using PIVlab and parallel computing. 
	This script creates a .mat file which contains : the computed velocity fields (u,v), PIV boxes position in pixels and tables of parameters (s and p) which 
	are used for the DIC algorithm

	This script calls function PIVlab_commandline_parallel_stereo described below

  --> PIVlab_commandline_parallel_stereo.m
	This function takes into account several arguments : 
	 - directory = file where .tiff images are saved 
	 - reference = reference image if we want to perform DIC using a single image as a reference frame, if reference = '', then images separated by time step 
	Dt are compared
	 - N = number of frames to process, if N = 0, the function process all images that are in the directory
	 - i0 = first image from which we start DIC
	 - Dt = time step, detailed above
	 - b, detailed above
	 - s = table of parameters used to perform PIV
	 - p = table of parameters for images pre-processing 

	This function calls piv_analysis.m, described below. 

	Then, we perform a loop over all images. Images are preprocessed by function PIVlabe_preproc which take as arguments an image and
	the table of parameters p. In our case we only use a contrast enhancement algorithm (CLAHE). Contrast is locally increased within a window of finite size,
	defined in table of parameters p ('CLAHE size').

	We then use piv_FFTmulti (described below) thanks to which we compute displacement field (u0,v0) between the two compared images. 

  --> piv_analysis.m
	Perform DIC algorithm on a pair of images (filename1, filename2). Processing is described in the Matlab function
	 



### Post-processing of raw displacement fields 
  --> Main_data_structuration.m
	This script order the raw data as well as the parameters used in PIVlabe into a matlab structure thanks to functions 
	PIV_banquise_postprocessing.m and genere_structure_banquise.m
	
	Different fields, parameters and time related to the drone flight are added to the structure which is saved under the format '*_total_processed.mat'

	From this first structure, a new structure, containing scaled velocity fields can be created using functuons 
	scaling_structure.m and scaling_structure_oblique.m. 
	These functions enable to scale the measured displacement fields. In the case of an oblique view, pitch angle not equal to 90°, vertical displacement field 
	is computed. This scaled version of the structure is saved under the name '*_scaled.mat'.

	Finally all parameters used during the process are saved in a .txt file thanks to the function saving_parameters.m



  --> PIV_banquise_postprocessing.m
	Post-process raw displacement fields obtained by function PIVlab_commandline_parallel_stereo.m
	This script takes as arguments : 
	- raw displacement fields (u,v)
	- parameters tables (s and p)
	- W, window size, previously described
	- N, number of frames to post-process, if N = 0, all frames must be post-processed

	It returns a filtered velocity field (u_filt,v_filt). 
	This script applies a median filter to the displacement fields and then replace outliers by interpolated values. Methods are described in 
	article 'Universal outlier detection for PIV data' Westerweel & Scarano (2005) https://doi.org/10.1007/s00348-005-0016-6


  --> genere_structure_banquise.m 
	Create a matlab structure (equivalent to a Python dictionnary) from a .mat file created by code automaticcycle_banquise_drone.m. 
	All inputs and outputs are defined in the function
	
	It is possible to supress PIV boxes on the image boundaries by using parameter a. This can be useful to always have all boxes of the same size. 


  --> scaling_structure.m / scaling_structure_oblique.m
	
	These function enable the scaling of measured velocity fields. 
	New fields are added to the pre-existing structure : 
		+ spatial coordinates x and y in meters, 
		+ experimental time t and international time t_UTC
		+ Vx, Vy and Vz, scaled velocity fields measured by PIV

	In the case of a vertical view (camera pitch angle of 90°), horizontal velocity fields Vx and Vy are scaled. In the case of an oblique view, 
	a vertical velocity field Vz is computed and scaled, horizontal fields are not to scale in this particular case. 












