#!/usr/bin/env python

#stdlib imports
import argparse
import os.path
import ConfigParser
import sys
import copy
import warnings

#third party imports
import numpy as np
from pyproj import Proj
from pagerio import shake
from pagerio import gmt,shapefile,esri
from shapely import geometry
from pagerutil import text

#set default grid resolution here - later add it to config files
DEFAULT_RES = 0.0083 #decimal degrees
    
def createGrids(global_grids,covshp,predictors,outfolder,ename,resolution=DEFAULT_RES):
    eventfolder = os.path.join(outfolder,ename)
    #figure out the minimum bounding box of the coverage data set
    shpobj = shapefile.PagerShapeFile(covshp)
    if not shpobj.hasIndex:
        shpobj.createShapeIndex()
    covbox = shpobj.bounds
    allbounds = [covbox]
    #figure out the minimum bounding box of the predictor grids
    for predictor_name,predictorvalue in predictors.iteritems():
        predictorfile,predictoratts = predictorvalue
        #what kind of a file are you?
        #We'll support more file types later, for now just handle shakemaps
        if not predictorfile.endswith('.xml'):
            f,e = os.path.splitext(predictorfile)
            raise Exception,'Unsupported file type for predictor: %s' % e
        
        shakemap = shake.ShakeGrid(predictorfile,variable=predictoratts[0].upper())
        shakebox = shakemap.getRange()
        allbounds.append(shakebox)
    xmin = -9999999999
    xmax = 9999999999
    ymin = xmin
    ymax = xmax
    for bounds in allbounds:
        txmin,txmax,tymin,tymax = bounds
        if txmin > xmin:
            xmin = txmin
        if txmax < xmax:
            xmax = txmax
        if tymin > ymin:
            ymin = tymin
        if tymax < ymax:
            ymax = tymax
    
    if xmin >= xmax or ymin >= ymax:
        raise Exception,'Your data sets do not overlap!'

    #make a reference geodict (clip everything else to this)
    ncols = int(round((xmax-xmin)/resolution))
    nrows = int(round((ymax-ymin)/resolution))
    #we need to ensure that the resulting grid is still equal to or inside the bounds
    #we just calculated.
    txmax = xmin + ncols*resolution
    while txmax > xmax:
        ncols -= 1
        txmax = xmin + ncols*resolution

    tymax = ymin + nrows*resolution
    while tymax > ymax:
        nrows -= 1
        tymax = ymin + nrows*resolution
        
    geodict = {'xmin':xmin,'xmax':xmax,'ymin':ymin,'ymax':ymax,
               'xdim':resolution,'ydim':resolution,
               'nrows':nrows,'ncols':ncols}    

    #create the output directory
    if not os.path.isdir(outfolder):
        os.makedirs(outfolder)
    
    #make a GMT grid from the shapefile data at the resolution specified - no clipping yet
    covgrid,ispoint = makeCoverageGrid(shpobj,resolution)
    covgrid.interpolateToGrid(geodict,method='nearest')
    covname = os.path.join(outfolder,'coverage.grd')
    print 'Saving coverage grid...'
    covgrid.save(covname)
    

    #clip out all of the global grids
    for gridkey,gridfilename in global_grids.iteritems():
        bounds = (xmin-(resolution*2),xmax+(resolution*2),ymin-(resolution*2),ymax+(resolution*2))
        if gridfilename.endswith('grd'):
            grid = gmt.GMTGrid(grdfile=gridfilename,bounds=bounds)
        else:
            #make sure we get a larger area than desired, so we can resample
            tmpgrid = esri.EsriGrid(gridfilename)
            tmpgrid.load(bounds=bounds)
            grid = gmt.GMTGrid()
            grid.loadFromGrid(tmpgrid)
        grid.interpolateToGrid(geodict)
        print 'Saving grid %s...' % gridkey
        grid.save(os.path.join(outfolder,gridkey+'.grd'))
        
    #Now save out the desired ShakeMap layers
    for predictor_name,predictorvalue in predictors.iteritems():
        predictorfile,predictoratts = predictorvalue
        #what kind of a file are you?
        #We'll support more file types later, for now just handle shakemaps
        if not predictorfile.endswith('.xml'):
            f,e = os.path.splitext(predictorfile)
            raise Exception,'Unsupported file type for predictor: %s' % e

        for variable in predictoratts:
            shakemap = shake.ShakeGrid(predictorfile,variable=variable.upper())
            shakemap.interpolateToGrid(geodict)
            gmtgrid = gmt.GMTGrid()
            gmtgrid.geodict = shakemap.geodict.copy()
            gmtgrid.griddata = shakemap.griddata.copy()
            print 'Saving grid %s...' % variable
            gmtgrid.save(os.path.join(outfolder,variable+'.grd'))
    
    return eventfolder

def makeCoverageGrid(covshp,resolution):
    ispoint = False
    if covshp.shapeType == 'polygon':
        covgrid = makePolygonGrid(covshp,resolution)
    elif covshp.shapeType == 'polygon':
        covgrid = makePointGrid(covshp,resolution)
        ispoint = True
    return (covgrid,ispoint)

def getBestUTM(cx):
    starts = np.arange(-180,180,6)
    zone = np.where((cx > starts) < 1)[0].min()
    return zone

def makePolygonGrid(covshp,resolution):
    covgrid = makeEmptyGrid(covshp.bounds,resolution,np.double)
    xmin,xmax,ymin,ymax = covshp.bounds
    utmzone = getBestUTM(xmin + (xmax-xmin)/2.0)
    utmproj = Proj(proj='utm',zone=utmzone,ellps='WGS84')
    shpcounter = 0
    progress_unit = text.roundToNearest(covshp.nShapes/10,1000)
    for shape in covshp.getShapes():
        if not shpcounter % progress_unit:
            print 'Processing shape %i of %i' % (shpcounter,covshp.nShapes)
        shapelon = shape['x']
        shapelat = shape['y']
        #sometimes these polygons are in pieces separated by nan values
        #in this case, we need to create a multipolygon object
        #do we have this kind of polygon?
        #if so, let's loop over each piece and do the cell area calculations
        if np.any(np.isnan(shapelat)):
            inan = np.isnan(shapelat).nonzero()[0].tolist()
            inan.append(len(shapelat))
            startidx = 0
            polylist = []
            for nanidx in inan:
                shapex,shapey = utmproj(shapelon[startidx:nanidx],shapelat[startidx:nanidx])
                slidepoly = geometry.Polygon(zip(shapex,shapey))
                countgrid(covgrid,shapelon,shapelat,slidepoly,resolution,utmproj)
                startidx = nanidx + 1
        #if not, just do calculations on simple polygon
        else:
            shapex,shapey = utmproj(shapelon,shapelat)
            slidepoly = geometry.Polygon(zip(shapex,shapey))
            countgrid(covgrid,shapelon,shapelat,slidepoly,resolution,utmproj)
       
        shpcounter += 1
    return covgrid

def countgrid(covgrid,shapelon,shapelat,slidepoly,resolution,utmproj):
    nrows,ncols = covgrid.griddata.shape
    sxmin = min(shapelon)
    sxmax = max(shapelon)
    symin = min(shapelat)
    symax = max(shapelat)
    uly,ulx = covgrid.getRowCol(symax,sxmin)
    lry,lrx = covgrid.getRowCol(symin,sxmax)
    if uly < 1:
        uly = 1
    if ulx < 1:
        ulx = 1
    if lry > nrows:
        lry = nrows
    if lrx > ncols:
        lrx = ncols
    for i in range(uly,lry+1):
        for j in range(ulx,lrx+1):
            clat,clon = covgrid.getLatLon(i,j)
            xmin = clon - resolution/2.0
            xmax = clon + resolution/2.0
            ymin = clat - resolution/2.0
            ymax = clat + resolution/2.0
            xcell = [xmin,xmin,xmax,xmax,xmin]
            ycell = [ymin,ymax,ymax,ymin,ymin]
            xcellutm,ycellutm = utmproj(xcell,ycell)
            cellpoly = geometry.Polygon(zip(xcellutm,ycellutm))
            try:
                if not cellpoly.intersects(slidepoly):
                    continue
            except:
                intpoly = cellpoly.intersection(slidepoly)
            cellfraction = slidepoly.area/cellpoly.area
            covgrid.griddata[i,j] += covgrid.griddata[i,j] + cellfraction

# def makePolygonGrid(covshp,resolution):
#     covgrid = makeEmptyGrid(covshp.bounds,resolution,np.double)
#     nrows,ncols = covgrid.griddata.shape
#     xmin,xmax,ymin,ymax = covshp.bounds
#     utmzone = getBestUTM(xmin + (xmax-xmin)/2.0)
#     utmproj = Proj(proj='utm',zone=utmzone,ellps='WGS84')
#     #this might be *really* slow...
#     for i in range(0,nrows):
#         print 'Row %i of %i...' % (i,nrows)
#         for j in range(0,ncols):
#             cxmin = xmin + j*resolution
#             cxmax = cxmin + resolution
#             cymax = ymax - i*resolution
#             cymin = cymax - resolution
#             lonbox = [cxmin,cxmax,cxmax,cxmin,cxmin]
#             latbox = [cymin,cymin,cymax,cymax,cymin]
#             #get all shapes intersecting cell
#             shapes = covshp.getShapesByBoundingBox((cxmin,cxmax,cymin,cymax))
#             if not len(shapes):
#                 continue
#             #union these shapes and project
#             unitedshape = projectUnitedShape(shapes,utmproj) #returns Shapely polygon
#             xbox,ybox = utmproj(lonbox,latbox)
#             cellshape = geometry.Polygon(zip(xbox,ybox))
#             #we only want the part of the united polygons that is inside this cell
#             unitedshape = unitedshape.intersection(cellshape)
#             fraction = unitedshape.area/cellshape.area
#             if fraction > 0:
#                 pass
#             covgrid.griddata[i,j] = fraction
#     return covgrid
            
def projectUnitedShape(shapes,utmproj):
    shape1lon = shapes[0]['x']
    shape1lat = shapes[0]['y']
    outx,outy = utmproj(shape1lon,shape1lat)
    outpoly = geometry.Polygon(zip(outx,outy))
    for i in range(1,len(shapes)):
        loncomp = shapes[i]['x']
        latcomp = shapes[i]['y']
        xcomp,ycomp = utmproj(loncomp,latcomp)
        polycomp = geometry.Polygon(zip(xcomp,ycomp))
        outpoly = outpoly.union(polycomp)
    return outpoly

def makePointGrid(covshp,resolution):
    covgrid = makeEmptyGrid(covshp.bounds,resolution,np.uint8)
    for feature in covshp.getShapes():
        x = feature['x']
        y = feature['y']
        row,col = covgrid.getRowCol(y,x)
        covgrid.griddata[row,col] = 1
    return covgrid

def makeEmptyGrid(covbox,resolution,dtype):
    covgrid = gmt.GMTGrid()
    ulx = covbox[0]
    uly = covbox[3]
    lrx = covbox[1]
    lry = covbox[2]
    ncols = int(np.ceil((lrx-ulx)/resolution))
    nrows = int(np.ceil((uly-lry)/resolution))
    xmin = ulx
    xmax = ulx + resolution*ncols
    ymax = uly
    ymin = uly - resolution*nrows
    geodict = {'xmin':ulx,'xmax':xmax,'ymin':ymin,'ymax':ymax,'xdim':resolution,'ydim':resolution}
    covgrid.geodict = geodict.copy()
    covgrid.griddata = np.zeros((nrows,ncols),dtype=dtype)
    return covgrid
    
def parseEvent(eventfile):
    reqsections = ['COVERAGE','PREDICTORS','ATTRIBUTES','OUTPUT']
    config = ConfigParser.RawConfigParser()
    config.read(eventfile)
    sections = config.sections()
    if not set(reqsections) <= set(sections): #is at least every required section present?
        raise Exception,'Incomplete config file - must have the following sections defined: %s' % str(reqoptions)
    missing = []
    covshp = config.get('COVERAGE','coverage')
    if not os.path.isfile(covshp):
        missing.append(covshp)
    predictors = dict(config.items('PREDICTORS'))
    for pgrid in predictors.values():
        if not os.path.isfile(pgrid):
            missing.append(pgrid)
    #make sure that the attribute names match the predictor names
    attnames = config.options('ATTRIBUTES')
    if not set(predictors.keys()) <= set(attnames):
        errstr = 'You must define a set of attributes for every predictor variable.'
        if len(missing):
            errstr += 'Also, the following files do not exist: %s' % str(missing)
        raise Exception,errstr
    for key,value in predictors.iteritems():
        attlist = config.get('ATTRIBUTES',key).split(',')
        predictors[key] = (value,attlist)

    if 'name' not in config.options('OUTPUT'):
        raise Exception,'Missing "name" field in OUTPUT section'
    ename = config.get('OUTPUT','name')
    return (covshp,predictors,ename)    

def writeConfig(global_config,outfolder,configfile):
    config = ConfigParser.RawConfigParser()
    config.add_section('GRIDS')
    for key,value in global_config.iter_items():
        config.set('GRIDS',key,value)
    config.add_section('OUTPUT')
    config.set('OUTPUT','folder',outfolder)
    p,f = os.path.split(configfile)
    if not os.path.isdir(p):
        os.makedirs(p)
    with open(configfile, 'wb') as f:
        config.write(f)

def readConfig(configfile):
    config = ConfigParser.RawConfigParser()
    config.read(configfile)
    global_grids = dict(config.items('GRIDS'))
    outfolder = config.get('OUTPUT','folder')
    return (global_grids,outfolder)

def configure(configfile):
    print 'You will now be prompted for the names and locations of the global grid files'
    print 'on your system.  When you are done entering these, hit enter when prompted for the'
    print 'name of the next grid.'
    global_config = {}
    while True:
        key = raw_input('Enter name of global grid (i.e., geology): ')
        if key.strip() == '':
            break
        value = raw_input('Enter path to global grid (i.e., /Users/user/geology.grd): ')
        if value.strip() == '':
            print "No grid path entered - let's try again"
            continue
        global_config[key] = value
    outfolder = raw_input('Enter the folder where all output should be written: ')
    writeConfig(global_config,outfolder,configfile)
    print 'Your config file has been written to %s' % configfile

def main(parser,args):
    configfile = os.path.join(os.path.expanduser('~'),'.lsprocess','lsprocess.cfg')
    hasconfig = os.path.isfile(configfile)
    if not hasconfig and not args.configure:
        print
        print 'No config file "%s" found.  \nRun with configure option to create.' % configfile
        print
        parser.print_help()
        sys.exit(1)
    if args.configure:
        if hasconfig:
            prompt = '''You already have a config file at %s.  
            Are you sure you want to re-configure? y/[n] ''' % configfile
            ans = raw_input(prompt)
            if ans.strip() == '' or ans.strip().lower() == 'n':
                sys.exit(0)
            configure(configfile)
        
    #read the global config
    global_grids,outfolder = readConfig(configfile) #returns a dictionary just like global_config above

    #validate/parse the event cfg file
    try:
        covshp,predictors,ename = parseEvent(args.eventfile)
    except Exception,msg:
        print 'There is something wrong with your event file.  See errors below.'
        print msg
        sys.exit(1)

    print 'Your coverage data set: %s' % covshp
    print 'Your global grids:'
    for key,value in global_grids.iteritems():
        print '\t%s = %s' % (key,value)
    print 'Your non-global predictor variables:'
    for key,value in predictors.iteritems():
        gridfile = value[0]
        gridattrs = ','.join(value[1])
        print '\t%s = %s (%s)' % (key,gridfile,gridattrs)
    print 'Your output will be written to:'
    eventfolder = os.path.join(outfolder,ename)
    print '\t%s' % eventfolder
    eventfolder = createGrids(global_grids,covshp,predictors,eventfolder,ename)
    
if __name__ == '__main__':
    #warnings.filterwarnings("ignore", "*matplotlib*")
    parser = argparse.ArgumentParser()
    parser.add_argument("eventfile", type=str,nargs='?',
                        help="A config file specifying event-specific input")
    parser.add_argument("-c", "--configure", action="store_true",
                        help="Create global config file")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="increase output verbosity")
    args = parser.parse_args()
    if args.configure:
        configfile = os.path.join(os.path.expanduser('~'),'.lsprocess','lsprocess.cfg')
        print 'Configuration tools not yet implemented.'
        print 'Create your global config file by hand here:\n'
        print configfile
        print
        sys.exit(0)
    main(parser,args)
        
    
    

