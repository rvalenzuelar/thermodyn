Thermodynamic functions included:

* saturation vapor pressure (Lowe 1980)
* saturation vapor pressure (Bohren and Albrecht 1998)
* saturation mixing ratio (Bohren and Albrecht 1998)
* dew point (Lawrence 2005)
* relative humidity (Lawrence 2005)
* virtual temperature (Saucier 1989)
* lifting condensation level temperature (Bolton 1980)
* potential temperature (Wallace and Hobbs 2006)
* potential temperature (Bolton 1980)
* equivalent potential temperature (Wallace and Hobbs 2006)
* equivalent potential temperature (Bolton 1980)
* Brunt-Vaisala frequency unsaturated (Bohren and Albrecht 1998)
* Brunt-Vaisala frequency saturated (Durran and Klemp 1982)


Need installed:

* numpy 
* pandas

Function arguments:

K = temperature in Kelvin 
C = temperature in deg Celsius
theta =  potential temperature [K]
Dewp = dew point temperature [C]
relh = relative humidity [%]
mixing_ratio = water vapor mixing ratio [kg/kg]
mb = pressure in milibar
hPa = pressure in hectopascal
agl_m = altitude above ground level (i.e. sounding altitude)
m = altitude above sea level (i.e. surface station)


Examples using ipython:

In [1]: import Thermodyn as thermo

In [2]: thermo.relative_humidity(C=20.0, Dewp=15.4)
Out[3]: array(77.0)

In [4]: thermo.relative_humidity(C=[20.0,25.0,23.1], Dewp=[15.4,25.0,23.2])
Out[5]: array([  77.,  100.,  100.])

In [6]: thermo.relative_humidity(K=[293.15,298.0,296.25], Dewp=np.array([288.55,299.0,299.0])-273.15)
Out[7]: array([  77,  100.  ,  100.  ])