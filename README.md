## About ACOLITE
ACOLITE combines the atmospheric correction algorithms for aquatic applications of Landsat and Sentinel-2 developed at RBINS.

It allows simple and fast processing of Landsat (5/7/8) and Sentinel-2 (A/B) images for coastal and inland water applications. Features include generation of RGB images before and after atmospheric correction, atmospheric correction of water bodies and extraction of rectangular regions of interest (defined by bounding coordinates). Level 2 outputs are surface reflectance (ρs=Rrs⋅π) and derived products that can be saved as PNG maps and geolocated datasets in a NCDF (NetCDF) file. The atmospheric correction is image based and needs no external inputs.

The algorithm was presented in Vanhellemont and Ruddick 2018, [Atmospheric correction of metre-scale optical satellite data for inland and coastal water applications](https://www.sciencedirect.com/science/article/pii/S0034425718303481) and Vanhellemont 2019, [Adaptation of the dark spectrum fitting atmospheric correction for aquatic applications of the Landsat and Sentinel-2 archives](https://doi.org/10.1016/j.rse.2019.03.010). As of October 2020, the default settings as suggested in Vanhellemont 2020, [Sensitivity analysis of the dark spectrum fitting atmospheric correction for metre- and decametre-scale satellite imagery using autonomous hyperspectral radiometry](https://doi.org/10.1364/OE.397456) are used, i.e. excluding the SWIR bands and using a new sky reflectance correction.


ACOLITE development was funded by the Belgian Science Policy Office STEREO program under contracts SR/37/135 (JELLYFOR project) and SR/00/325 (PONDER project), and by the European Community's Seventh Framework Programme (FP7/2007-2013) under grant agreement n° 606797 (HIGHROC project).

**ACOLITE is provided by RBINS as an experimental tool, without explicit or implied warranty. Use of the program is at your own discretion and risk.**

## Distribution
ACOLITE is distributed as a binary package on the [REMSEM page](http://odnature.naturalsciences.be/remsem/software-and-data/acolite) and supported on the [ACOLITE forum](http://odnature.naturalsciences.be/remsem/acolite-forum/).

This git repository serves as a distribution of the source code and is aimed only at experienced users.

## Dependencies
ACOLITE is coded in Python 3, and requires the following Python packages to run with all functionality:`matplotlib scipy pyproj gdal netcdf4 pyhdf requests statsmodels basemap pillow scikit-image pyresample`

## Installation
* cd into a suitable directory and clone the git repository: `git clone https://github.com/acolite/acolite`
* cd into the new acolite directory `cd acolite`
* run `python launch_acolite.py`
