## About ACOLITE
ACOLITE combines the atmospheric correction algorithms for aquatic applications of Landsat and Sentinel-2 developed at RBINS.

It allows simple and fast processing of Landsat (5/7/8) and Sentinel-2 (A/B) images for coastal and inland water applications. Features include generation of RGB images before and after atmospheric correction, atmospheric correction of water bodies and extraction of rectangular regions of interest (defined by bounding coordinates). Level 2 outputs are surface reflectance (ρs=Rrs⋅π) and derived products that can be saved as PNG maps and geolocated datasets in a NCDF (NetCDF) file. The atmospheric correction is image based and needs no external inputs. A full publication describing the new dark spectrum fitting algorithm is forthcoming.

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
