# Wind

## Requirements

### Python packages

You will need the following packages:
```
xarray
numpy
matplotlib
```
You can use a pip env or a conda env.

### Files

You need to download from the server `canots` the files.
Copying one by one can be tedious, so one way to speed up this is to gather all
files of interest in a directory on the server, then copy every files from this
directory onto the local machine.

If you want to get all data from the Trisonica that was mounted
on the mobile mast, just type (zsh)
```
ls **/Mat_portatif/Trisonica/*.txt
```
or use the little script (bsh)
```
./gather_files.sh
```
or you can use rsync (not tested)
```
rsync -av user@hostname:'/path/to/root/**/Mat_portatif/Trisonica/*.txt' /home/jacqhugo/BicWin26/terrain/data_trisonica_portable/
```


## Reading anemometers

### Principle of measure

An anemometer measure the wind speed in the three directions of space, and the
temperature at a high sampling rate (from 5Hz to 20Hz). A sound wave is emitted
to a receptor at some distance. From delay between the speed of sound (computed
separately) and the speed of the sound wave a wind speed is deduced. A large (20cm) 
distance between emitter and receiver and a high sampling rate (20Hz) allow for
a precise detection of the turbulence. The closer to the ground, the smaller and
higher frequency you need to adapt to the size of the coherent structures
(scales as the distance to the wall).

### Code example

Take a look at `first_look.py` for a example of data processing for the wind
sensors (Trisonica and Thies). You can write the path to the file and name to
the file inside the script. Then, simply run
```
python first_look.py
```

The script mainly uses the *parse_anemo_to_xarray* function inside
`read_anemo.py`. We use one function for both anemometers as the global
structure of the data file is the same. The difference is in the line validation
function (*validate_line*) and the actual line parsing function (*read_data*)
that both require the `kind` argument.

### Eddy covariance

To do

### Monin-Obukov (log law)

To do

## Reading IR120

### Principle of measure

We measure a flux in long wave radiations. The sensor (thermopile) is emitting
heat, influencing its own measurement. To counter this, a resistor (the
thermistor) with a know thermal behavior is set up to be able to deduce how much
the flux measure is modified by the body of the sensor. Alongside the flux, we
need to measure the resistance of the thermistor from which we deduce the
temperature of the thermistor. Then, we can also correct for the IR flux
received by the sensor from the ambiant air. Finally, with a known emissivity we
get the surface temperature.

The surface temperature is used as:
- context of the day: clouds, coldness of the day, . ..
- reference temperature for the computation of turbulent fluxes


### Code example

You can look at `test_marmotte.py` where you can build a temperature evolution
plot from a data file. This script is build upon the `read_marmotte.py` script
that parses the data from the file. To convert the raw data to a temperature,
the script `sensorIR120.py` is used, with calibration constant from the
manufacturer (A,B,C,X,Y,Z).

```
python test_marmotte.py
```

## Method

A detailed method is in the code files

## References

> Stull, R. B. (1988). An introduction to boundary layer meteorology. Springer Dordrecht.  

> Kaimal, J. C., Finnigan, J. J., & Kaimal, J. C. (1994). Atmospheric boundary layer flows: Their structure and measurement. Oxford University Press.  

> Lee, X., Massman, W., & Law, B. (Eds.). (2005). Handbook of Micrometeorology (Vol. 29). Springer Netherlands. https://doi.org/10.1007/1-4020-2265-4  
Foken, T., & Mauder, M. (2024). Micrometeorology (Third edition). Springer.  

> Aubinet, M., Vesala, T., & Papale, D. (Eds.). (2012). Eddy Covariance: A Practical Guide to Measurement and Data Analysis. Springer Netherlands. https://doi.org/10.1007/978-94-007-2351-1  

> Sicart, J. E. (2021). Contribution à l’étude des flux d’énergie en surface de glaciers de montagne dans les Andes tropicales et dans les Alpes. https://hal.science/tel-03337891v1  


