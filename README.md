TODO
____
 - Automatically projecting coverage
 - Handling coverage formats other than shapefile
 - Handling (projecting as well) any other predictor data in vector format

Introduction
------------

lsprocess is a project to convert data in various forms into a package
of raster data files suitable for use in logistic regression for the
PAGER Secondary Hazards Project.  

Data can be broken up into three
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

Usage:
---------------------
A script called "secsample" will be installed in your path.  This script should be called from the command line with one argument,
the name of a config file whose format is described below.  The script will create a training and (possibly) a testing CSV file
in the designated output directory.

Sample INI file
--------
<pre>
[SAMPLING]
coverage = /path/to/hazard_coverage/shapefile #required, must be in decimal degrees
dx = 10.0 #sampling resolution in meters, required
cb = 0.5  #forced hazard/no-hazard class balance, optional.  Number specifies the fraction of hazard pixels to sample
nmax = 1e9 #optional maximum number of possible yes/no sample points (usually set to avoid memory issues)
nsamp = 1e5 #optional number of total hazard and no-hazard sample points to collect.
testpercent = 0.5 #Fraction of sampled points to be used for testing (1-testpercent) fraction will be used for training. Optional, defaults to 0
extent = xmin,xmax,ymin,ymax OR convex #geographic extent within which to sample data.  Four numbers are interpreted as bounding box, the word convex will be interpreted to mean a convex hull.  Default (not specified) will mean the bounding box of the hazard coverage.

h1 = 100.0 #Minimum buffer size for sampling non-hazard points when input coverage takes the form of points.
h2 = 300.0 #Maximum buffer size for sampling non-hazard points when input coverage takes the form of points.

[PREDICTORS]
layername = /path/to/predictor_grid/or/predictor_shapefile #inputs can be ESRI or GMT format grids, or shapefiles.  Must be in decimal degrees.
layername_sampling = nearest #optional grid sampling method (nearest or linear will be supported)
layername_attribute = attribute #required for shapefiles - the attribute of each shape to choose as sample.

shakemap = /path/to/grid.xml #required parameters specifying path to ShakeMap input.  All ground motion values (mmi,pga,pgv,psa03,psa10,psa30) will be sampled.
shakemap_uncertainty = /path/to/uncertainty.xml #optional path to ShakeMap uncertainty grid.  All error columns corresponding to ground motions will be sampled.

[OUTPUT]
folder = /path/to/desired/output/location #required, where all data frames, output plots, etc. will be written
basename = eqname #base name to assign to all output files (eqname_testing.dat, eqname_training.dat, etc.)
</pre>

