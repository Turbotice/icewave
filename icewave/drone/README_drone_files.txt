# Details about files that can be extracted from Drone DJI Mavic 3 Pro

## Different ypes of files

- Videos, format .MP4
- Videos aformat .LRF (Low Resolution File), can be converted to .MP4 to get a low resolution file
- fichier format .SRT (SubRip Subtitle) file, used to add subtitles to videos
It is possible to extract from these files the Following parameters for each frame of the video performed : 
	- experimental time and UTC time at which the frame has been captured 
	- camera propreties : (iso, shutter, focal length)
	- drone GPS location (1e-5 precision) and relative altitude in meter (1e-1 meter precision), relative to the home point altitude
- photo au format .JPG 
- flightrecords.txt : encrypted binary txt file… It is hard to decode this document. The only way to get these data is by using PhantomHelp LogViewer (https://www.phantomhelp.com/LogViewer/upload/) or the designed
software : flightreader (see section below)

## How to organize files ? 

1. Sort flights by date, exemple : 2023/0130/ 
2. Set the experiment index number of the flight : 01,  02,  03 ...
3. Gather all different flies for a single experiment : MP4, LRF, SRT, ...
4. Use convert_multivideo2tiff.py to generate '.tiff' images from a video
5. The video is ready to be analyzed (PIV, Optical flow, Stereo imaging)


## How to get flightrecords.txt data ?

1. Connect the controler to a desktop thanks to a USB cable
2. Gather flight files : " flightrecords" for each flight performed
3. Use the website PhantomHelp (https://www.phantomhelp.com/LogViewer/upload/), or use the software flightreader, load a flightrecord.txt
4. Download the generated .csv and .kml files 
5. Gather these files in the corresponding folder 

More info on this forum : https://phantompilots.com/threads/flight-log-retrieval-and-analysis-guide.147496/


## Relevant parameters recorded by the drone 

From .SRT files :
1. Local time (set on UAV) of each frame
2. Latitude and longitude for each frame
3. Relative altitude of the drone (home point altitude chosen as a reference) for each frame (in meters)

From Flightrecords csv files : All parameters listed below are recorded every 0.1 second
1. Local time
2. Latitude and longitude
3. OSD.height (ft) UAV height (relative to home point) in feets
4. OSD.altitude(ft) UAV altitude = home point altitude + OSD.height (in feets)
5. OSD.yaw [360], UAV orientation with respect to North direction (0 is North), computed by the IMU (in degrees)
6. GIMABL.pitch, Camera angle to the horizontal (in degrees)
7. GIMBAL.yaw [360], camera orientation with respect to North direction (0 is North), computed by IMU (in degrees)
8. CAMERA.isPhoto , boolean whether if the camera is capturing pictures or not. Not very reliant..
9. CAMERA.isVideo, boolean whether if the camera is capturing videos or not. Not very reliant..

Many other parameters are recorded and stored in these csv files  


