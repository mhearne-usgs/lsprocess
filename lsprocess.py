#!/usr/bin/env python

#stdlib import
import sys
import os.path
import argparse
import configparser

#third party
import numpy as np
import fiona
import fiona.crs as crs
from neicio.shake import ShakeGrid
from neicio.esri import EsriGrid
from neicio.gmt import GMTGrid
from neicutil.matutil import repmat
import rasterio
from affine import Affine
from rasterio import features
from scipy.io import netcdf

def getFileResolution(gridfile):
    xdim = None
    ydim = None
    if gridfile.endswith('grd'):
        cdf = netcdf.netcdf_file(gridfile)
        xvarname = None
        if 'x' in list(cdf.variables.keys()):
            xvarname = 'x'
            yvarname = 'y'
        elif 'lon' in list(cdf.variables.keys()):
            xvarname = 'lon'
            yvarname = 'lat'
        if xvarname is not None: #at least two forms of COARDS-compliant netcdf files...
            xvar = cdf.variables[xvarname].data.copy()
            yvar = cdf.variables[yvarname].data.copy()

            #do some QA on the x and y data
            dx = np.diff(xvar)
            dy = np.diff(yvar)

            isXConsistent = np.abs(1 - np.max(dx)/np.min(dx)) < 0.01
            isYConsistent = np.abs(1 - np.max(dx)/np.min(dx)) < 0.01
            if not isXConsistent or not isYConsistent:
                raise Exception('X or Y cell dimensions are not consistent!')

            #assign x/y resolution
            xdim = np.mean(dx)
            ydim = np.mean(dy)
        cdf.close()
    else:
        gridbase,gridext = os.path.splitext(gridfile)
        headerfile = gridbase + '.hdr'
        f = open(headerfile,'rt')
        for line in f.readlines():
            parts = line.split()
            key = parts[0]
            if key == 'xdim':
                xdim = float(parts[1])
            if key == 'ydim':
                ydim = float(parts[1])
            if key == 'cellsize':
                xdim = float(parts[1])
                ydim = float(parts[1])
                break
        f.close()
    return (xdim,ydim)
        

def sampleGrid(gridfile,geodict):
    xdim,ydim = getFileResolution(gridfile)
    resolution = np.mean([xdim,ydim])
    xmin = geodict['xmin']
    xmax = geodict['xmax']
    ymin = geodict['ymin']
    ymax = geodict['ymax']
    #resolution = geodict['xdim']
    bounds = (xmin-(resolution*2),xmax+(resolution*2),ymin-(resolution*2),ymax+(resolution*2))
    if gridfile.endswith('grd'):
        grid = GMTGrid(grdfile=gridfile,bounds=bounds)
    else:
        tmpgrid = EsriGrid(gridfile)
        tmpgrid.load(bounds=bounds)
        grid = GMTGrid()
        grid.loadFromGrid(tmpgrid)
    grid.interpolateToGrid(geodict)
    return grid

def getShape(bbox,resolution):
    xmin = bbox[0]
    ymax = bbox[3]
    xmax = bbox[1]
    ymin = bbox[2]
    ncols = int(np.ceil((xmax-xmin)/resolution))
    nrows = int(np.ceil((ymax-ymin)/resolution))
    xmax = xmin + ncols*resolution
    ymin = ymax - nrows*resolution
    outbbox = [xmin,xmax,ymin,ymax]
    return (nrows,ncols,outbbox)

def makeCoverageGrid(covshp,geodict):
    shapes = fiona.open(covshp)
    geoms = []
    for shape in shapes:
        geoms.append(shape['geometry'])
    shapes.close()
    outshape = (geodict['nrows'],geodict['ncols'])
    transform = Affine.from_gdal(geodict['xmin'],geodict['xdim'],0.0,geodict['ymax'],0.0,-geodict['ydim'])
    img = features.rasterize(geoms,out_shape=outshape,fill=0,
                             transform=transform,all_touched=True,
                             default_value=1)
    covgrid = GMTGrid()
    covgrid.geodict = geodict
    covgrid.griddata = np.int8(img.copy())
    return covgrid

def parseEvent(eventfile):
    reqsections = ['COVERAGE','PREDICTORS','ATTRIBUTES','OUTPUT']
    config = configparser.RawConfigParser()
    config.read(eventfile)
    sections = config.sections()
    if not set(reqsections) <= set(sections): #is at least every required section present?
        raise Exception('Incomplete config file - must have the following sections defined: %s' % str(reqoptions))
    missing = []
    covdict = {}
    fname = config.get('COVERAGE','coverage')
    covdict['filename'] = fname
    covpath,covfile = os.path.split(fname)
    covname,covext = os.path.splitext(covfile)
    prjfile = os.path.join(covpath,covname,'.prj')
    hasProjFile = os.path.isfile(prjfile)
    hasProjStr = config.has_option('COVERAGE','projstr')
    # if not hasProjFile and not hasProjStr:
    #     raise Exception,'You must have a .prj file or a proj4 string set in config file.'
    
    # covdict['filename'] = fname
    # if hasProjStr:
    #     covdict['projstr'] = config.get('COVERAGE','projstr')
    # else:
    #     source = fiona.open(fname)
    #     covdict['projstr'] = covdictcrs.to_string(source.crs)
    #     source.close()

    #get the format and bbox fields (if found)
    if config.has_option('COVERAGE','format'):
        covdict['format'] = config.get('COVERAGE','format')
    if config.has_option('COVERAGE','bbox'):
        covdict['bbox'] = [float(b) for b in config.get('COVERAGE','bbox').split()]

    #look for resolution field
    if config.has_option('COVERAGE','resolution'):
        covdict['resolution'] = float(config.get('COVERAGE','resolution'))
        
    predictors = dict(config.items('PREDICTORS'))
    for pgrid in list(predictors.values()):
        if pgrid.find('_projstr') > -1 or pgrid.find('_format') > -1:
            continue
        if not os.path.isfile(pgrid):
            missing.append(pgrid)
    #make sure that the attribute names match the predictor names
    attnames = config.options('ATTRIBUTES')
    if not set(predictors.keys()) <= set(attnames):
        errstr = 'You must define a set of attributes for every predictor variable.'
        if len(missing):
            errstr += 'Also, the following files do not exist: %s' % str(missing)
        raise Exception(errstr)
    for key,value in predictors.items():
        attlist = config.get('ATTRIBUTES',key).split(',')
        predictors[key] = (value,attlist)

    if 'name' not in config.options('OUTPUT'):
        raise Exception('Missing "name" field in OUTPUT section')
    ename = config.get('OUTPUT','name')
    return (covdict,predictors,ename)    

def readConfig(configfile):
    config = configparser.RawConfigParser()
    config.read(configfile)
    global_grids = dict(config.items('GRIDS'))
    outfolder = config.get('OUTPUT','folder')
    return (global_grids,outfolder)

def main(args):
    #read in global config file
    configfile = os.path.join(os.path.expanduser('~'),'.lsprocess','lsprocess.cfg')
    hasconfig = os.path.isfile(configfile)
    if not hasconfig:
        print()
        print('No config file "%s" found.' % configfile)
        print()
        sys.exit(1)
    global_grids,outfolder = readConfig(configfile) #returns a dictionary just like global_config above
    
    
    #read in event specific grid file
    try:
        covdict,predictors,ename = parseEvent(args.eventfile)
    except Exception as msg:
        print('There is something wrong with your event file.  See errors below.')
        print(msg)
        sys.exit(1)
    
    #construct output folder from global/event configs
    outfolder = os.path.join(outfolder,ename)
    if not os.path.isdir(outfolder):
        os.mkdir(outfolder)
    
    #look for bounding box and resolution in event config file, or get from shakemap
    bbox = None
    shakemap = ShakeGrid(predictors['shakemap'][0],'MMI')
    if 'bbox' in covdict:
        bbox = covdict['bbox']
    else:
        #bbox = shakemap.getRange()
        #default to the bounding box of the coverage data
        with fiona.open(covdict['filename']) as src:
            tbbox = src.bounds
            bbox = (tbbox[0],tbbox[2],tbbox[1],tbbox[3])
            
    if 'resolution' in covdict:
        resolution = covdict['resolution']
    else:
        resolution = shakemap.getGeoDict()['xdim']
    
    #get input coverage projection from event config OR from .prj file
    #projstr = covdict['projstr']
    
    #get format of coverage, check against list of supported fiona formats, read in data
    #we'll do other support later
    
    #if necessary, project coverage into lat/lon
    #skip projection for now as well

    #determine what the grid shape and (potentially) new bbox is given bbox and resolution
    nrows,ncols,bbox = getShape(bbox,resolution)
    #if the coverage dataset is larger than the ShakeMap, we need to make sure our output grid
    #is contained by the shakemap for interpolation purposes.
    shakebounds = shakemap.getRange()
    shakexdim,shakeydim = (shakemap.geodict['xdim'],shakemap.geodict['ydim'])
    xmin = max(bbox[0],shakebounds[0]+shakexdim*2)
    xmax = min(bbox[1],shakebounds[1]-shakexdim*2)
    ymin = max(bbox[2],shakebounds[2]+shakeydim*2)
    ymax = min(bbox[3],shakebounds[3]-shakeydim*2)
    geodict = {'xdim':resolution,'ydim':resolution,
               'xmin':xmin,'xmax':xmax,
               'ymin':ymin,'ymax':ymax,
               'nrows':nrows,'ncols':ncols}
    
    #rasterize projected coverage defined bounding box and resolution
    shpfile = covdict['filename']
    print('Creating coverage grid...')
    covgrid = makeCoverageGrid(shpfile,geodict)
    outgridfile = os.path.join(outfolder,'coverage.grd')
    print('Saving coverage to %s...' % outgridfile)
    covgrid.save(outgridfile)

    #make a grid of lat,lon values
    row = np.arange(0,nrows)
    col = np.arange(0,ncols)
    rows = repmat(row,ncols,1).T
    cols = repmat(col,nrows,1)
    lat,lon = covgrid.getLatLon(rows,cols)

    #create a list of arrays that we'll dump out to a text file when done
    vardict = {}
    vardict['coverage'] = covgrid.griddata.flatten()
    vardict['lat'] = lat.flatten()
    vardict['lon'] = lon.flatten()
        
    #subset shakemap and global grids using defined bounding box and resolution
    shakefile = predictors['shakemap'][0]
    variables = predictors['shakemap'][1]
    for var in variables:
        shakemap = ShakeGrid(shakefile,var.upper())
        shakemap.interpolateToGrid(geodict)
        gmtshake = GMTGrid()
        gmtshake.geodict = shakemap.geodict
        gmtshake.griddata = shakemap.griddata
        outshakefile = os.path.join(outfolder,'%s.grd' % var)
        print('Saving %s to %s...' % (var,outshakefile))
        gmtshake.save(outshakefile)
        vardict[var] = gmtshake.griddata.flatten()
        
    #write netcdf versions of coverage, shakemap, and global grids to output folder
    for gridname,gridfile in global_grids.items():
        if not os.path.isfile(gridfile):
            pass
        try:
            grid = sampleGrid(gridfile,geodict)
        except Exception as e:
            print('There was an error while sampling the "%s" grid "%s". - "%s"' % (gridname,gridfile,str(e)))
            
        outgridfile = os.path.join(outfolder,gridname+'.grd')
        print('Saving %s to %s...' % (gridname,outgridfile))
        grid.save(outgridfile)
        vardict[gridname] = grid.griddata.flatten()
        
    #create text file with columns of data for all predictor variables
    firstcols = ['lat','lon','coverage']
    outmat = np.zeros((nrows*ncols,len(vardict)))
    for i in range(0,len(firstcols)):
        col = firstcols[i]
        outmat[:,i] = vardict[col]
    colidx = i+1
    colnames = []
    for col,column in vardict.items():
        if col in firstcols:
            continue
        outmat[:,colidx] = vardict[col]
        colnames.append(col)
        colidx += 1

    colnames = firstcols + colnames
    m,n = outmat.shape
    datfile = os.path.join(outfolder,'%s.dat' % ename)
    print('Saving all variables to data file %s...' % datfile)
    f = open(datfile,'wt')
    f.write(','.join(colnames)+'\n')
    for i in range(0,m):
        line = ','.join('%.4f' % col for col in outmat[i,:])
        f.write(line+'\n')
    f.close()
    
    
if __name__ == '__main__':
    aparser = argparse.ArgumentParser()
    aparser.add_argument("eventfile", type=str,nargs='?',
                        help="A config file specifying event-specific input")
    pargs = aparser.parse_args()
    main(pargs)
    
    
    
    
    
    
