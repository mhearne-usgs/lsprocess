TODO
____
 - Automatically projecting coverage
 - Handling coverage formats other than shapefile
 - Handling (projecting as well) any other predictor data in vector format

Introduction
------------

lsprocess is a project to convert data in various forms into a package
of raster data files suitable for use in logistic regression for the
PAGER Secondary Hazards Project.  Data can be broken up into three
categories:

 * Global data - Global raster data sets which will be applied to all
   input events.
 * Coverage data - An event specific landslide coverage vector data
   set, in ShapeFile or XY (Lon/Lat) format.
 * Predictor data - Event specific raster/vector data sets which are
   used as predictor variables in logistic regression.

The lsprocess application will read a global config file, which will
contain the paths to the above referenced global raster data sets.  It
will then take as input a second *event-specific* config file which
will specify the coverage and predictor data (see above).

Installation and Dependencies
-----------------------------

This package depends on:
 * numpy, the fundamental package for scientific computing with Python. <a href="http://www.numpy.org/">http://www.numpy.org/</a>  
 * fiona, a Python library for reading/writing various vector data formats (shapefile)
 * rasterio, a Python library used in this case for rasterizing vector data.
 * neicio, a Python library for reading/writing various grid data formats (ShakeMap and GMT).
 * neicutil, a catch-all utility Python library (used in this case for the repmat function).

The best way to install numpy is to use one of the Python distributions described here:

<a href="http://www.scipy.org/install.html">http://www.scipy.org/install.html</a>

Anaconda and Enthought distributions have been successfully tested with smtools.

Most of those distributions should include <em>pip</em>, a command line tool for installing and 
managing Python packages.  You will use pip to install the other dependencies and smtools itself.  
 
You may need to open a new terminal window to ensure that the newly installed versions of python and pip
are in your path.

To install fiona:

pip install fiona

To install rasterio, follow the instructions found here:

<a href="https://github.com/mapbox/rasterio">https://github.com/mapbox/rasterio</a>

To install neicio:

pip install git+git://github.com/usgs/neicio.git

To install neicutil:

pip install git+git://github.com/usgs/neicutil.git

To install this package:

pip install git+git://github.com/mhearne-usgs/lsprocess.git

Uninstalling and Updating
-------------------------

To uninstall:

pip uninstall lsprocess

To update:

pip install -U git+git://github.com/mhearne-usgs/lsprocess.git

Sample GLOBAL INI file
--------
<pre>
[GRIDS]
#these can be any file format supported by GDAL (?)
geology = /Users/user/data/geology.flt
slopemin = /Users/user/data/slopemin.flt
slopemax = /Users/user/data/slopemax.flt
cti = /Users/user/data/cti.grd

[OUTPUT]
#each data run will be saved in a folder under here
folder = /Users/user/lsprocess/
</pre>

Sample event-specific INI file
--------
<pre>
[COVERAGE]
coverage = /Users/user/events/event1/coverage.shp
resolution = 0.008337 #defaults to the resolution of the ShakeMap

[PREDICTORS]
shakemap = /Users/user/events/event1/grid.xml #ShakeMap grid XML format

#It is possible to have shapefile vector data sets to use as predictor variables
geology_event = /Users/user/events/event1/geology.shp 

[ATTRIBUTES]

#here we're specifying the attributes of the predictor data sets that
#should be output (shakemap) or rasterized (vector)

shakemap = pga,pgv #list of layers to output

geology_event = rock_type #name of the shapefile attribute to rasterize

[OUTPUT]
name = event1
</pre>

Usage
--------

usage: lsprocess.py [-h] [eventfile]

positional arguments:
  eventfile   A config file specifying event-specific input

optional arguments:
  -h, --help  show this help message and exit
