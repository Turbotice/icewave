from PIL import Image
from PIL.ExifTags import TAGS

import matplotlib.pyplot as plt
import numpy as np

global folder
folder = '/Users/stephane/Documents/BicWin2026/Data_Exemples/IR_0211/'
# path to the image or video
imagename = folder + "DJI_20260212054323_0043_T.jpg"


# read the image data using PIL
im = plt.imread(imagename)
image = Image.open(imagename)


#plt.imshow(im)
#plt.show()

# extract EXIF data
exifdata = image.getexif()

# iterating over all EXIF data fields
for tag_id in exifdata:
    # get the tag name, instead of human unreadable tag id
    tag = TAGS.get(tag_id, tag_id)
    data = exifdata.get(tag_id)#.decode("utf-16")
    print(f"{tag:25}: {data}")  
