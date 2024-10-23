


global date
date = "0223'

############### Drones #############################
import icewave.field.drone as drone
#generate records dictionnary from srt and jpg files
#(jpg files are currently empty)
drone.get_records(date)

#Convert flightrecords in csv format to pkl format
drone.convert_flightrecords(date)

#find all mp4 files
#generate a .tiff image for the first frame of each movie
drone.get_mp4files(date,save=True)

############### Phones ##############################
import icewave.phone.analyse as analyse

# step1 : extract all data in .zip format
#analyse.process(date,1)

# step2 : find the measurement interval, generate a dataset in pkl format
analyse.process(date,2)

# step3 : compute average quantities, generate a summary csv file
analyse?process(date,3)


############### Multi instruments ###################
import icewave.field.multi_instruments as multi

#save the records file for all instruments
multi.save_records(date)
