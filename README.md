Introduction
------------

lsprocess is a project to convert data in various forms into a package of raster data files suitable
for use in logistic regression for the PAGER Secondary Hazards Project.  Data can be broken up into three categories:
 * Global data - Global raster data sets which will be applied to all input events.
 * Coverage data - An event specific landslide coverage vector data set, in ShapeFile or XY (Lon/Lat) format.
 * Predictor data - Event specific raster/vector data sets which are used as predictor variables in logistic regression.

The lsprocess application will read a global config file, which will contain the paths to the above referenced global 
raster data sets.  It will then take as input a second *event-specific* config file which will specify the coverage 
and predictor data (see above).

Sample GLOBAL INI file
--------
[GRIDS]
geology = /Users/user/data/geology.flt
slopemin = /Users/user/data/slopemin.flt
slopemax = /Users/user/data/slopemax.flt
cti = /Users/user/data/cti.grd

Sample event-specific INI file
--------
[COVERAGE]
coverage = /Users/user/events/event1/coverage.shp

[PREDICTORS]
shakemap = /Users/user/events/event1/grid.xml #ShakeMap grid XML format

#we can possibly have shapefile vector data sets to use as predictor variables
geology_event = /Users/user/events/event1/geology.shp 

[ATTRIBUTES]
#here we're specifying the attributes of the predictor data sets that should be output (shakemap) or rasterized (vector)
shakemap = pga,pgv #list of layers to output
geology_event = rock_type #name of the shapefile attribute to rasterize

[OUTPUT]
name = event1

Dummy Examples
--------

    >>> foo = bar()




