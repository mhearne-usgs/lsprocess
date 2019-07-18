TODO
____
 - Handling custom data extents
 - Handling linear grid sampling method

Introduction
------------

lsprocess is a project to convert data in various forms into a package
of raster data files suitable for use in logistic regression for the
PAGER Secondary Hazards Project.  

Installation and Dependencies
-----------------------------

Currently, to install the code and its dependencies requires use of git.

Do the following (If you already have Anaconda or Miniconda installed, skip steps #2-3):

  1. git clone https://github.com/mhearne-usgs/lsprocess
  2. cd lsprocess
  3. ./install.sh
  4. conda activate pysample


Usage:
---------------------
A script called "secsample" will be installed in your path.  This script should be called from the command line with one argument,
the name of a config file whose format is described below.  The script will create a training and (possibly) a testing CSV file
in the designated output directory.

It is also possible to use secsample to get information about the sampling that will be done using the config file.  Below is a sample 
of the type of output you would get:

<pre>
secsample -c earthquake.ini
At a resolution of 30.0 meters, the input shapefile /Users/mhearne/data/landslide/northridge/northridge_dd.shp would have:
	2,482 rows
	2,769 columns
	6,872,658 total possible samples
	A class balance of 0.39% hazard pixels
	Estimated number of hazard pixels: 26,474
	Estimated upper bound for nsamp: 6,872,540
</pre>

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

#Optional Boolean (0 or 1) indicating whether presence of polygon in each grid cell is enough to turn that
#into a yes pixel.  Setting this to 0 presumes that the dx is relatively large, such
#that creating a grid at that resolution will not tax the resources of the system.
#Set this to 1 if doing high-resolution sampling of data, 0 if coarse (i.e., grid cells of 1km).
touch_center = 0 

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

