#!/usr/bin/env python

#stdlib imports
import os.path
from functools import partial
from collections import OrderedDict

import numpy as np
import fiona
from shapely.geometry import Polygon,shape,MultiPoint,Point,mapping
import matplotlib.pyplot as plt
import pyproj
import pandas as pd
from shapely.ops import transform
from mapio.grid2d import Grid2D
from mapio.gmt import GMTGrid
from mapio.gdal import GDALGrid
from mapio.shake import ShakeGrid
from mapio.geodict import GeoDict
from rasterio.transform import Affine
import rasterio

def getPointsInCircum(r,n=100,h=0,k=0):
    #h = x coordinate of center
    #k = y coordinate of center
    #http://stackoverflow.com/questions/8487893/generate-all-the-points-on-the-circumference-of-a-circle
    points = [(np.cos(2*np.pi/n*x)*r,np.sin(2*np.pi/n*x)*r) for x in xrange(0,n+1)]
    x,y = zip(*points)
    x = np.array(x)
    y = np.array(y)
    x += h
    y += k
    return (x,y)

def createCirclePolygon(h,k,r,dx):
    D = 10.0
    theta = 2 * np.arccos((r-(dx/D))/r)
    npoints = int(360.0/theta)
    x,y = getPointsInCircum(r,n=npoints,h=h,k=k)
    p = Polygon(zip(x,y))
    return p

def affine_from_corner(ulx, uly, dx, dy):
    return Affine.translation(ulx, uly)*Affine.scale(dx, -dy)

def getNoSampleGrid(yespoints,xvar,yvar,dx,h1,h2):
    '''Return the grid from which "no" pixels can successfully be sampled.
    :param yespoints:
      Sequence of (x,y) points (meters) where landslide/liquefaction was observed.
    :param xvar:
      Numpy array of centers of columns of sampling grid.
    :param yvar:
      Numpy array of centers of rows of sampling grid.
    :param dx:
      Sampling resolution in x and y (meters).
    :param h1:
      Minimum buffer size for sampling non-hazard points.
    :param h2:
      Maximum buffer size for sampling non-hazard points.
    :returns:
      Grid of shape (len(yvar),len(xvar)) where 1's represent pixels from which 
      "no" values can be sampled.
    '''
    shp = (len(xvar),len(yvar))
    west = xvar.min() - dx/2.0 #??
    north = yvar.max() + dx/2.0 #??
    affine = affine_from_corner(west,north,dx,dx)
    donuts = []
    holes = []
    for h,k in yespoints:
        donut = createCirclePolygon(h,k,h2,dx)
        hole = createCirclePolygon(h,k,h1,dx)
        donuts.append(donut)
        holes.append(hole)
    donutburn = ((mapping(g), 1) for g in donuts)
    holeburn = ((mapping(g), 2) for g in holes)
    #we only want those pixels set where the polygon encloses the center point
    alltouched = False 
    donutimg = rasterio.features.rasterize(
                donutburn,
                out_shape=shp,
                transform=affine,
                all_touched=alltouched)
    holeimg = rasterio.features.rasterize(
                holeburn,
                out_shape=shp,
                transform=affine,
                all_touched=alltouched)
    holeimg[holeimg == 0] = 1
    holeimg[holeimg == 2] = 0
    sampleimg = np.bitwise_and(donutimg,holeimg)
    return sampleimg

def getFileType(filename):
    """
    Determine whether input file is a shapefile or a grid (ESRI or GMT).
    :param filename:
      String path to candidate filename.
    :returns:
      String, one of 'shapefile','grid','unknown'.
    """
    fname,fext = os.path.splitext(filename)
    dbf = fname + '.dbf'
    ftype = 'unknown'
    if os.path.isfile(dbf):
        ftype = 'shapefile'
    else:
        try:
            fdict = GMTGrid.getFileGeoDict(filename)
            ftype = 'grid'
        except Exception,error:
            try:
                fdict,xvar,yvar = GDALGrid.getFileGeoDict(filename)
                ftype = 'grid'
            except:
                pass
    return ftype

def getProjectedShapes(shapes,xmin,xmax,ymin,ymax):
    """
    Take a sequence of geographic shapes and project them to a bounds-centered orthographic projection.
    :param shapes:
      Sequence of shapes, as read in by fiona.collection()
    :param xmin:
      Eastern boundary of all shapes
    :param xmax:
      Western boundary of all shapes
    :param ymin:
      Southern boundary of all shapes
    :param ymax:
      Northern boundary of all shapes
    :returns:
       - Input sequence of shapes, projected to orthographic
       - PyProj projection object used to transform input shapes
    """
    latmiddle = ymin + (ymax-ymin)/2.0
    lonmiddle = xmin + (xmax-xmin)/2.0
    projstr = '+proj=ortho +datum=WGS84 +lat_0=%.4f +lon_0=%.4f +x_0=0.0 +y_0=0.0' % (latmiddle,lonmiddle)
    proj = pyproj.Proj(projparams=projstr)
    project = partial(
        pyproj.transform,
        pyproj.Proj(proj='latlong', datum='WGS84'),
        proj)

    pshapes = []
    for tshape in shapes:
        if tshape['geometry']['type'] == 'Polygon':
            pshapegeo = shape(tshape['geometry'])
        else:
            pshapegeo = shape(tshape['geometry'])
        pshape = transform(project, pshapegeo)
        pshapes.append(pshape) #assuming here that these are simple polygons

    return (pshapes,proj)

def getYesPoints(pshapes,proj,dx,nmax):
    """
    Collect x/y coordinates of all points within hazard coverage polygons at desired resolution.
    :param pshapes:
      Sequence of orthographically projected shapes.
    :param proj:
      PyProj projection object used to transform input shapes
    :param dx:
      Float resolution of grid at which to sample points
    :param nmax:
      Threshold maximum number of points in total data mesh.
    :returns:
      - numpy 2-D array of X/Y coordinates inside hazard polygons.
      - number of rows of resulting mesh
      - number of columns of resulting mesh
      - numpy array of x coordinate centers of columns
      - numpy array of y coordinate centers of rows
      - 1D array of indices where yes pixels are located (use np.unravel_index to unpack to 2D array)
    """
    mxmin = 9e10
    mxmax = -9e10
    mymin = 9e10
    mymax = -9e10
    for pshape in pshapes:
        pxmin,pymin,pxmax,pymax = pshape.bounds
        if pxmin < mxmin:
            mxmin = pxmin
        if pxmax > mxmax:
            mxmax = pxmax
        if pymin < mymin:
            mymin = pymin
        if pymax > mymax:
            mymax = pymax

    xvar = np.arange(mxmin,mxmax+dx,dx)
    yvar = np.arange(mymin,mymax+dx,dx)
    ncols = len(xvar)
    nrows = len(yvar)
    if nmax is not None:
        if ncols*nrows > nmax:
            aspect = ncols/nrows
            ncols = np.sqrt(nmax*aspect)
            nrows = nmax/ncols
            ncols = int(ncols)
            nrows = int(nrows)
            #re-calculate dx here...
            tdx = (mxmax-mxmin)/ncols
            tdy = (mymax-mymin)/nrows
            dx = np.max(tdx,tdy)
            xvar = np.arange(mxmin,mxmax+dx,dx)
            yvar = np.arange(mymin,mymax+dx,dx)

    #Get the "yes" points to sample from
    yespoints = []
    idx = []
    shapeidx = 0
    if pshapes[0].type == 'Polygon':
        #loop over shapes, projecting each one, then get the sample points
        for pshape in pshapes:
            if not shapeidx % 1000:
                print 'Searching polygon %i of %i' % (shapeidx,len(pshapes))
            shapeidx += 1
            pxmin,pymin,pxmax,pymax = pshape.bounds
            leftcol = np.where((pxmin - xvar) >= 0)[0].argmax()
            rightcol = np.where((xvar - pxmax) >= 0)[0][0]
            bottomrow = np.where((pymin - yvar) >= 0)[0].argmax()
            toprow = np.where((yvar - pymax) >= 0)[0][0]
            xp = np.arange(xvar[leftcol],xvar[rightcol]+dx,dx)
            yp = np.arange(yvar[bottomrow],yvar[toprow]+dx,dx)
            xmesh,ymesh = np.meshgrid(xp,yp)
            xy = zip(xmesh.flatten(),ymesh.flatten())
            for point in xy:
                ix = np.where(xvar == point[0])[0][0]
                iy = np.where(yvar == point[1])[0][0]
                if pshape.contains(Point(point)):
                    yespoints.append(point)
                    idx.append(np.ravel_multi_index((iy,ix),(nrows,ncols),mode='raise',order='C'))
    else:
        yespoints = []
        for pshape in pshapes:
            yespoints.append(pshape.coords[0])
            
    return (np.array(yespoints),nrows,ncols,xvar,yvar,idx)

def sampleFromFile(shapefile,dx=10.0,nmax=None,testPercent=0.0,classBalance=None,extent=None,Nsamp=100,h1=100.0,h2=300.0):
    """
    Sample yes/no test and training pixels from shapefile input.
    :param shapefile:
       path to shapefile, presumed to be decimal degrees.
    :param dx:
       resolution of sampling in X and Y (meters)
    :param nmax:
      if not None, maximum allowed number of mesh points in X and Y together (nrows*ncols).  Overrides dx.
    :param testPercent:
      Fraction of total samples to put into output testing data frame. (training frame gets 1-testpercent).  Default is 0.0.
    :param classBalance:
      If None, uses class balance of data (raises exception for points).  If specified, is the fraction of "yes" pixels, where "yes"+"no" = 1
    :param extent:
      If none, use the bounding box of the data in the shapefile. 
    :param Nsamp:
      Number of total (yes+no) sample points.
    :returns:
      - sequence of XY coordinates for:
         - YesTestPoints
         - YesTrainPoints
         - NoTestPoints
         - NoTrainPoints
      - numpy array of mesh column centers
      - numpy array of mesh row centers
      - PyProj object defining orthographic projection of xy points
    """
    #read the shapes in from the file
    f = fiona.collection(shapefile,'r')
    shapes = list(f)
    bounds = f.bounds
    
    f.close()

    return sampleFromShapes(shapes,bounds,dx=dx,nmax=nmax,testPercent=testPercent,
                            classBalance=classBalance,extent=extent,Nsamp=Nsamp,
                            h1=h1,h2=h2)

def sampleYes(array,N):
    """
    Sample without replacement N points from an array of XY coordinates.
    :param array:
      2D numpy array of XY points
    :param N:
      int number of points to sample without replacement from input array.
    :returns:
      Tuple of (sampled points, unsampled points)
    """
    #array is a Mx2 array of X,Y points
    m,n = array.shape
    allidx = np.arange(0,m)
    sampleidx = np.random.choice(allidx,size=N,replace=False)
    nosampleidx = np.setxor1d(allidx,sampleidx)
    sampled = array[sampleidx,:]
    notsampled = array[nosampleidx,:]
    return (sampled,notsampled)

def sampleNo(xvar,yvar,N,avoididx):
    """
    Sample from pixels in mesh, excluding yes pixels and already sampled no pixels.
    :param xvar:
      Numpy array of centers of all columns in mesh.
    :param yvar:
      Numpy array of centers of all rows in mesh.
    :param N:
      Number of no pixels to sample.
    :param avoididx:
      1D array of indices from mesh that should NOT be sampled from.  Initially this will be the array
      of indices where the yes pixels are.
    :returns:
      Randomly chosen list of tuples of (x,y) coordinate points that are outside polygons.
    """
    allidx = np.arange(0,len(xvar)*len(yvar)) #flattened array of all indices in mesh
    noidx = np.setxor1d(allidx,avoididx) #allidx - avoididx
    #noidx = np.array(list(set(allidx) - set(avoididx)))
    nosampleidx = np.random.choice(noidx,size=N,replace=False)
    newavoididx = np.sort(np.hstack((avoididx,nosampleidx)))
    rowidx,colidx = np.unravel_index(nosampleidx,(len(yvar),len(xvar)))
    samples = []
    for row,col in zip(rowidx,colidx):
        xp = xvar[col]
        yp = yvar[row]
        samples.append((xp,yp))

    return (samples,newavoididx)

def sampleFromShapes(shapes,bounds,dx=10.0,nmax=None,testPercent=1.0,classBalance=None,extent=None,Nsamp=100,h1=100.0,h2=300.0):
    """
    Sample yes/no test and training pixels from shapefile input.
    :param shapes:
       Sequence of projected shapes.
    :param dx:
       resolution of sampling in X and Y (meters)
    :param nmax:
      if not None, maximum allowed number of mesh points in X and Y together (nrows*ncols).  Overrides dx.
    :param testPercent:
      Fraction of total samples to put into output testing data frame. (training frame gets 1-testpercent).  Default is 0.0.
    :param classBalance:
      If None, uses class balance of data (raises exception for points).  If specified, is the fraction of "yes" pixels, where "yes"+"no" = 1
    :param extent:
      If none, use the bounding box of the data in the shapefile. 
    :param Nsamp:
      Number of total (yes+no) sample points.
    :returns:
      - sequence of XY coordinates for:
         - YesTestPoints
         - YesTrainPoints
         - NoTestPoints
         - NoTrainPoints
      - numpy array of mesh column centers
      - numpy array of mesh row centers
      - PyProj object defining orthographic projection of xy points
    """
    xmin,ymin,xmax,ymax = bounds
    shptype = shapes[0]['geometry']['type']
    if shptype not in ['Point','Polygon']:
        raise Exception('Only polygon and point data types supported!')

    #Get the shapes projected into an orthographic projection centered on the data 
    pshapes,proj = getProjectedShapes(shapes,xmin,xmax,ymin,ymax)

    if shptype != 'Polygon':
        if classBalance is None:
            raise Exception('Class balance *must* be selected when input data are points.')
    else:
        #what is the class balance (assuming polygons)
        if classBalance is None:
            classBalance = getClassBalance(pshapes,bounds,proj)

    #get the "yes" sample points
    yespoints,nrows,ncols,xvar,yvar,yesidx = getYesPoints(pshapes,proj,dx,nmax)

    #Calculations of how many training and test points are the same for points and polygons.
    #Also sampling of yes points is the same regardless of vector type
    Nmesh = nrows*ncols

    #Nsamp may not have been set - until we support custom extents by polygon, just assume that default Nsamp is = Nmesh
    if Nsamp is None:
        print 'Assuming total number of samples is the same as the number of pixels in sampling grid.  This may not work...'
        Nsamp = Nmesh
    
    NyesTot = len(yespoints)
    NnoTot = Nmesh - NyesTot
    NyesSampTest = int(Nsamp * classBalance * testPercent)
    NyesSampTrain = int(Nsamp * classBalance * (1-testPercent))
    YesSampTot = NyesSampTest + NyesSampTrain
    ratio = NyesTot/float(YesSampTot)
    if YesSampTot > NyesTot:
        raise Exception('Your total number of desired "yes" sample pixels is greater than the number available.')
    NnoSampTest = int(Nsamp*(1-classBalance)*testPercent)
    NnoSampTrain = int(Nsamp*(1-classBalance)*(1-testPercent))
    NoSampTot = NnoSampTest + NnoSampTrain
    if NoSampTot > NnoTot:
        raise Exception('Your total number of desired "no" sample pixels is greater than the number available.')
    YesTestPoints,RemainingYesPoints = sampleYes(yespoints,NyesSampTest)
    YesTrainPoints,RemainingYesPoints = sampleYes(RemainingYesPoints,NyesSampTrain)

    #Sampling of "no" points differs between points and polygons
    if shptype == 'Point':
        #for point data, create a boolean grid located on xvar/yvar
        #create a donut around each pixel containing a yes point(s), where 1 is assigned to
        #pixels more than h1 meters from yes point and less than h2 meters from yes point.
        nosampleimg = getNoSampleGrid(yespoints,xvar,yvar,dx,h1,h2)
        NoTestPoints,nosampleimg,sampleidx = sampleNoPoints(nosampleimg,NnoSampTest,xvar,yvar)
        NoTrainPoints,nosampleimg,sampleidx = sampleNoPoints(nosampleimg,NnoSampTrain,xvar,yvar)
    else:
        if extent is None: #we're using the default bounding box of the coverage data
            NoTestPoints,nosampleidx = sampleNo(xvar,yvar,NnoSampTest,yesidx)
            NoTrainPoints,nosampleidx = sampleNo(xvar,yvar,NnoSampTrain,nosampleidx)
        else:
            raise Exception('Custom extents not yet supported')  

    #project all of the point data sets back to lat/lon
    YesTestPoints = projectBack(YesTestPoints,proj)
    YesTrainPoints = projectBack(YesTrainPoints,proj)
    NoTestPoints = projectBack(NoTestPoints,proj)
    NoTrainPoints = projectBack(NoTrainPoints,proj)
    return (YesTestPoints,YesTrainPoints,NoTestPoints,NoTrainPoints,xvar,yvar,pshapes,proj)

def sampleNoPoints(sampleimg,N,xvar,yvar):
    '''Sample from our "no" sample grid, where that grid contains 1s where no pixels should be sampled from, and 0s where they should not.
    :param sampleimg:
      Grid of shape (len(yvar),len(xvar)) where 1's represent pixels from which 
      "no" values can be sampled.
    :param N:
      Number of pixels to sample (without replacement from sampleimg.
    :param xvar:
      Numpy array of centers of columns of sampling grid.
    :param yvar:
      Numpy array of centers of rows of sampling grid.
    :returns:
      - nopoints sequence of x,y tuples representing "no" samples.
      - Modified version of input sampleimg, with nopoints pixels set to 0.
    '''
    #get N points from sampleimg without replacement
    #avoid nosampleidx indices
    #return an sequence of X,Y tuples from those indices
    npixels = len(xvar)*len(yvar)
    allidx = np.arange(0,npixels)
    sampleimg = sampleimg.flatten() #flatten out the input image
    sampleidx = np.random.choice(allidx[sampleimg==1],size=N,replace=False)
    sampleidx.sort()
    sampleimg[sampleidx] = 0
    sampleimg.shape = (len(yvar),len(xvar))
    sampley,samplex = np.unravel_index(sampleidx,sampleimg.shape)
    xp = xvar[samplex]
    yp = yvar[sampley]
    nopoints = zip(xp,yp)
    return (nopoints,sampleimg,sampleidx)
    
    

def projectBack(points,proj):
    """
    Project a 2D array of XY points from orthographic projection to decimal degrees.
    :param points:
      2D numpy array of XY points in orthographic projection.
    :param proj:
      PyProj object defining projection.
    :returns:
      2D numpy array of Lon/Lat coordinates.
    """
    mpoints = MultiPoint(points)
    project = partial(
        pyproj.transform,
        proj,
        pyproj.Proj(proj='latlong', datum='WGS84'))
    gmpoints = transform(project, mpoints)
    coords = []
    for point in gmpoints.geoms:
        x,y = point.coords[0]
        coords.append((x,y))
    coords = np.array(coords)
    return coords

def plotPoints(shapes,YesTestPoints,YesTrainPoints,NoTestPoints,NoTrainPoints,filename):
    """
    Plot yes/no sample points and polygons.
    :param shapes:
      shapes (decimal degrees) as read in by fiona.collection()
    :param YesTestPoints:
      numpy 2D array of testing hazard sample points, decimal degrees.
    :param YesTrainPoints:
      numpy 2D array of training hazard sample points, decimal degrees. 
    :param NoTestPoints:
      numpy 2D array of testing non-hazard sample points, decimal degrees.
    :param NoTrainPoints:
      numpy 2D array of training non-hazard sample points, decimal degrees. 
    """
    #plot the "yes" sample points and the projected polygons
    figure = plt.figure(figsize=(8,8))
    plt.hold(True)
    for shape in shapes:
        px,py = zip(*shape['geometry']['coordinates'])
        px = list(px)
        py = list(py)
        px.append(px[0])
        py.append(py[0])
        plt.plot(px,py,'b')
    yestestx,yestesty = zip(*YesTestPoints)
    yestrainx,yestrainy = zip(*YesTrainPoints)
    plt.plot(yestestx,yestesty,'r.')
    plt.plot(yestrainx,yestrainy,'g.')
    notestx,notesty = zip(*NoTestPoints)
    notrainx,notrainy = zip(*NoTrainPoints)
    plt.plot(notestx,notesty,'cx')
    plt.plot(notrainx,notrainy,'mx')
    plt.legend(['Polygon 1 boundaries','Polygon 2 boundaries','Yes Test','Yes Train','No Test','No Train'],numpoints=1)
    plt.title('Sample points')
    plt.savefig(filename)

def getClassBalance(pshapes,bounds,proj):
    """
    Get native class balance of projected shapes, assuming a rectangular bounding box.
    :param pshapes:
      Sequence of projected shapely shapes
    :param bounds:
      Desired bounding box, in decimal degrees.
    :param proj:
      PyProj object defining orthographic projection of shapes.
    :returns:
      Float fraction of hazard polygons (area of hazard polygons/total area of bbox) 
    """
    xmin,ymin,xmax,ymax = bounds
    bpoly = Polygon([(xmin,ymax),
                     (xmax,ymax),
                     (xmax,ymin),
                     (xmin,ymin)])
    project = partial(
        pyproj.transform,
        pyproj.Proj(proj='latlong', datum='WGS84'),
        proj)
    bpolyproj = transform(project,bpoly)
    totalarea = bpolyproj.area
    polyarea = 0
    for pshape in pshapes:
        polyarea += pshape.area

    return polyarea/totalarea

def sampleShapeFile(shapefile,xypoints,attribute):
    """
    Open a shapefile (decimal degrees) and get the attribute value at each of the input XY points. Slower than sampling grids.
    :param shapefile:
      ESRI shapefile (decimal degrees) of predictor variable.
    :param xypoints:
      2D numpy array of XY points, in decimal degrees.
    :param attribute:
      String name of attribute to sample in each of the shapes.
    :returns:
      1D array of attribute values at each of XY points.
    """
    xmin = np.min(xypoints[:,0])
    xmax = np.max(xypoints[:,0])
    ymin = np.min(xypoints[:,1])
    ymax = np.max(xypoints[:,1])
    #xypoints should be projected back to lat/lon
    f = fiona.collection(shapefile,'r')
    tshapes = list(f.items(bbox=(xmin,ymin,xmax,ymax)))
    shapes = []
    for fid,shape in tshapes:
        shapes.append(shape)
    f.close()
    return sampleShapes(shapes,xypoints,attribute)

def subsetShapes(shapefile,bounds):
    """
    Return the subset of shapes from a shapefile within a given bounding box. (Can be slow).
    :param shapefile:
      ESRI shapefile (decimal degrees) of predictor variable.
    :param bounds:
      Bounding box (decimal degrees) to use for subset (xmin,ymin,xmax,ymax)
    """
    xmin,ymin,xmax,ymax = bounds
    #xypoints should be projected back to lat/lon
    f = fiona.collection(shapefile,'r')
    tshapes = list(f.items(bbox=(xmin,ymin,xmax,ymax)))
    shapes = []
    for fid,shape in tshapes:
        shapes.append(shape)
    f.close()
    return shapes

def sampleShapes(shapes,xypoints,attribute):
    """
    Get the attribute value at each of the input XY points for sequence of input shapes (decimal degrees). Slower than sampling grids.
    :param shapes:
      sequence of shapes (decimal degrees) of predictor variable.
    :param xypoints:
      2D numpy array of XY points, in decimal degrees.
    :param attribute:
      String name of attribute to sample in each of the shapes.
    :returns:
      1D array of attribute values at each of XY points.
    """
    samples = []
    for x,y in xypoints:
        for tshape in shapes:
            polygon = shape(tshape['geometry'])
            point = Point(x,y)
            if polygon.contains(point):
                sample = tshape['properties'][attribute]
                samples.append(sample)
    return np.array(samples)

def sampleGridFile(gridfile,xypoints,method='nearest'):
    """
    Sample grid file (ESRI or GMT format) at each of a set of XY (decimal degrees) points.
    :param gridfile:
      Name of ESRI or GMT grid format file from which to sample values.
    :param xypoints:
      2D numpy array of XY points, decimal degrees.
    :param method:
      Interpolation method, either 'nearest' or 'linear'.
    :returns:
      1D numpy array of grid values at each of input XY points.
    """
    xmin = np.min(xypoints[:,0])
    xmax = np.max(xypoints[:,0])
    ymin = np.min(xypoints[:,1])
    ymax = np.max(xypoints[:,1])
    gridtype = None
    try:
        fdict = GMTGrid.getFileGeoDict(gridfile)
        gridtype = 'gmt'
    except Exception,error:
        try:
            fdict,xvar,yvar = GDALGrid.getFileGeoDict(gridfile)
            gridtype = 'esri'
        except:
            pass
    if gridtype is None:
        raise Exception('File "%s" does not appear to be either a GMT grid or an ESRI grid.' % gridfile)
    xmin = xmin - fdict['dx']*3
    xmax = xmax + fdict['dx']*3
    ymin = ymin - fdict['dy']*3
    ymax = ymax + fdict['dy']*3
    bounds = (xmin,xmax,ymin,ymax)
    if gridtype == 'gmt':
        fgeodict = GMTGrid.getFileGeoDict(gridfile)
    else:
        fgeodict,xvar,yvar = GDALGrid.getFileGeoDict(gridfile)
    dx,dy = (fgeodict['dx'],fgeodict['dy'])
    try:
        sdict = GeoDict({'xmin':bounds[0],'xmax':bounds[1],
                        'ymin':bounds[2],'ymax':bounds[3],
                        'dx':dx,'dy':dy,
                        'ny':2,'nx':2},preserve='dims')
    except:
        pass
    if gridtype == 'gmt':
        grid = GMTGrid.load(gridfile,samplegeodict=sdict,resample=False,method=method,doPadding=True)
    else:
        grid = GDALGrid.load(gridfile,samplegeodict=sdict,resample=False,method=method,doPadding=True)

    return sampleFromGrid(grid,xypoints)

def sampleFromGrid(grid,xypoints,method='nearest'):
    """
    Sample 2D grid object at each of a set of XY (decimal degrees) points.
    :param grid:
      Grid2D object at which to sample data.
    :param xypoints:
      2D numpy array of XY points, decimal degrees.
    :param method:
      Interpolation method, either 'nearest' or 'linear'.
    :returns:
      1D numpy array of grid values at each of input XY points.
    """
    samples = []
    for lon,lat in xypoints:
        sample = grid.getValue(lat,lon,method=method)
        samples.append(sample)

    return np.array(samples)

def sampleFromShakeMap(shakefile,layer,xypoints):
    """
    Sample ShakeMap grid file at each of a set of XY (decimal degrees) points.
    :param shakefile:
      Grid2D object at which to sample data.
    :param xypoints:
      2D numpy array of XY points, decimal degrees.
    :returns:
      1D numpy array of grid values at each of input XY points.
    """
    shakegrid = ShakeGrid.load(shakefile,fixFileGeoDict='corner')
    return sampleFromMultiGrid(shakegrid,layer,points)

def sampleFromMultiGrid(multigrid,layer,xypoints):
    """
    Sample MultiGrid object (like a ShakeGrid) at each of a set of XY (decimal degrees) points.
    :param multigrid:
      MultiGrid object at which to sample data.
    :param xypoints:
      2D numpy array of XY points, decimal degrees.
    :returns:
      1D numpy array of grid values at each of input XY points.
    """
    if layer not in multigrid.getLayerNames():
        raise Exception('Layer %s not found in grid' % layer)
    hazgrid = multigrid.getLayer(layer)
    return sampleFromGrid(hazgrid,xypoints)

def getDataFrames(sampleparams,shakeparams,predictors,outparams):
    """
    Return Pandas training and testing data frames containing sampled data from hazard coverage, ShakeMap, and predictor data sets.
    :param sampleparams:
      Dictionary with at least these values:
        - coverage: Name of hazard coverage shapefile (decimal degrees). Required.
        - dx: Float desired sample resolution, and can be overridden by nmax, below (meters).  Required.
        - cb: Desired class balance, i.e., fraction of sampled points that should be from hazard polygons. Optional for polygons, Required for points.
        - nmax: Maximum number of possible yes/no sample points (usually set to avoid memory issues). Optional.
        - nsamp: Number of total hazard and no-hazard sample points to collect.  Required.
        - testpercent: Fraction of sampled points to be used for testing (1-testpercent) will be used for training. Optional, defaults to 0
        - extent: xmin,xmax,ymin,ymax OR convex #geographic extent within which to sample data.  Four numbers are interpreted as bounding box, the word convex will be interpreted to mean a convex hull.  Default (not specified) will mean the bounding box of the hazard coverage. Optional.
        - h1: Minimum buffer size for sampling non-hazard points when input coverage takes the form of points. Optional for polygons, required for points.
        - h2: Maximum buffer size for sampling non-hazard points when input coverage takes the form of points. Optional for polygons, required for points.
    :param shakeparams:
      Dictionary with at least these values:
        - shakemap: Name of shakemap file to use for sampling hazard values. Required.
        - shakemap_uncertainty: Name of shakemap uncertainty file to use for sampling hazard uncertainty values. Optional.
    :param predictors:
      Dictionary with at least these values:
        - layername: Path to ESRI shapefile, or grid in GMT or ESRI format which represents predictor data. Required.
        - layername_sampling: 'nearest' or 'linear', optional for grids, not used for shapefiles.
        - layername_attribute: Name of attribute in shapefile which should be sampled at hazard/non-hazard points.  Required for points.
    :param outparams:
      Dictionary with at least these values:
        - folder: Name of folder where all output (data frames, plots) will be written.  Will be created if does not exist. Required.
        - basename: The name that will be included in all output file names (i.e., northridge_train.csv). Required.
    :returns:
      Tuple of (training,testing) Pandas data frames. 
    """
    coverage = sampleparams['coverage']
    f = fiona.collection(coverage,'r')
    cbounds = f.bounds
    f.close()
    dx = sampleparams['dx']
    cb = sampleparams['cb']
    nmax = sampleparams['nmax']
    nsamp = sampleparams['nsamp']
    testpercent = sampleparams['testpercent']
    extent = sampleparams['extent']
    h1 = sampleparams['h1']
    h2 = sampleparams['h2']

    yestest,yestrain,notest,notrain,xvar,yvar,pshapes,proj = sampleFromFile(coverage,dx=dx,nmax=nmax,testPercent=testpercent,
                                                                            classBalance=cb,extent=extent,Nsamp=nsamp,h1=h1,h2=h2)

    traincolumns = OrderedDict()
    testcolumns = OrderedDict()

    if (100-testpercent) > 0:
        traincolumns['lat'] = np.concatenate((yestrain[:,1],notrain[:,1]))
        traincolumns['lon'] = np.concatenate((yestrain[:,0],notrain[:,0]))
        traincolumns['coverage'] = np.concatenate((np.ones_like(yestrain[:,1]),np.zeros_like(notrain[:,1])))

    if testpercent > 0:
        testcolumns['lat'] = np.concatenate((yestest[:,1],notest[:,1]))
        testcolumns['lon'] = np.concatenate((yestest[:,0],notest[:,0]))
        testcolumns['coverage'] = np.concatenate((np.ones_like(yestest[:,1]),np.zeros_like(notest[:,1])))
    
    
    for predname,predfile in predictors.iteritems():
        ftype = getFileType(predfile)
        if ftype == 'shapefile':
            attribute = predictors[predname+'_attribute']
            shapes = subsetShapes(predfile,cbounds)
            yes_test_samples = sampleShapes(shapes,yestest,attribute)
            no_test_samples = sampleShapes(shapes,notest,attribute)
            yes_train_samples = sampleShapes(shapes,yestrain,attribute)
            no_train_samples = sampleShapes(shapes,notrain,attribute)
            testcolumns[predname] = np.concatenate((yes_test_samples,no_test_samples))
            traincolumns[predname] = np.concatenate((yes_train_samples,no_train_samples))
        elif ftype == 'grid':
            method = 'nearest'
            if predictors.has_key(predname+'_sampling'):
                method = predictors[predname+'_sampling']

            if testpercent > 0:
                yes_test_samples = sampleGridFile(predfile,yestest,method=method)
                no_test_samples = sampleGridFile(predfile,notest,method=method)
                testcolumns[predname] = np.concatenate((yes_test_samples,no_test_samples))

            if (100-testpercent) > 0:
                yes_train_samples = sampleGridFile(predfile,yestrain,method=method)
                no_train_samples = sampleGridFile(predfile,notrain,method=method)
                traincolumns[predname] = np.concatenate((yes_train_samples,no_train_samples))
        else:
            continue #attribute or sampling method key

    #sample the shakemap
    layers = ['mmi','pga','pgv','psa03','psa10','psa30']
    shakegrid = ShakeGrid.load(shakeparams['shakemap'],fixFileGeoDict='corner')
    for layer in layers:
        yes_test_samples = sampleFromMultiGrid(shakegrid,layer,yestest)
        no_test_samples = sampleFromMultiGrid(shakegrid,layer,notest)
        yes_train_samples = sampleFromMultiGrid(shakegrid,layer,yestrain)
        no_train_samples = sampleFromMultiGrid(shakegrid,layer,notrain)
        testcolumns[layer] = np.concatenate((yes_test_samples,no_test_samples))
        traincolumns[layer] = np.concatenate((yes_train_samples,no_train_samples))
        
    dftest = pd.DataFrame(testcolumns)
    dftrain = pd.DataFrame(traincolumns)

    return (dftrain,dftest)
    

def test_sampleArray():
    array = np.random.rand(10,2)
    #get 5 random non-replacement
    allrows = np.arange(0,10)
    sampleidx = np.random.choice(allrows,size=5,replace=False)
    print sampleidx

def _test_sample():
    shapes = []
    poly1 = [(34.25,34.5),
             (34.3,34.4),
             (34.4,34.4),
             (34.45,34.5),
             (34.5,34.7),
             (34.4,34.7),
             (34.3,34.65),
             (34.25,34.6)]
    poly2 = [(34.75,34.65),
             (34.75,34.35),
             (34.8,34.35),
             (34.85,34.3),
             (34.95,34.3),
             (34.95,34.6)]
    shp1 = {'id':0,
            'type':'Feature',
            'properties':OrderedDict([('value',1)]),
            'geometry':{'type':'Polygon','coordinates':poly1}}
    shp2 = {'id':1,
            'type':'Feature',
            'properties':OrderedDict([('value',2)]),
            'geometry':{'type':'Polygon','coordinates':poly2}}
    shapes.append(shp1)
    shapes.append(shp2)
    allverts = np.vstack((np.array(poly1),np.array(poly2)))
    xmin = allverts[:,0].min()
    xmax = allverts[:,0].max()
    ymin = allverts[:,1].min()
    ymax = allverts[:,1].max()
    bounds = (xmin,ymin,xmax,ymax)
    dx = 5500.0
    YesTestPoints,YesTrainPoints,NoTestPoints,NoTrainPoints,xvar,yvar,pshapes,proj = sampleFromShapes(shapes,bounds,dx=dx,Nsamp=50,testPercent=0.5)
    plotPoints(shapes,YesTestPoints,YesTrainPoints,NoTestPoints,NoTrainPoints,'output.png')

    #make up a geology vector data set
    geopoly1 = [(34.0,34.6),
                (34.0,34.0),
                (34.1,34.0),
                (34.15,34.1),
                (34.3,34.25),
                (34.35,34.4),
                (34.25,34.5),
                (34.25,34.6)]
    geopoly2 = [(34.25,34.6),
                (34.25,34.5),
                (34.35,34.4),
                (34.6,34.4),
                (34.6,34.55),
                (34.6,34.6)]
    geopoly3 = [(34.1,34.0),
                (34.85,34.0),
                (34.6,34.4),
                (34.35,34.4),
                (34.3,34.25),
                (34.15,34.1)]
    geopoly4 = [(34.85,34.0),
                (35.3,34.0),
                (35.3,34.6),
                (34.6,34.6),
                (34.6,34.4),
                (34.85,34.15)]
    value = 4
    idx = 0
    shapes = []
    for poly in [geopoly1,geopoly2,geopoly3,geopoly4]:
        shp = {'id':idx,
               'type':'Feature',
               'properties':OrderedDict([('value',value)]),
               'geometry':{'type':'Polygon','coordinates':poly}}
        shapes.append(shp)
        value += 1
        idx += 1
    
    #make up a grid data set
    geodict = {'xmin':33.5,'xmax':35.5,'ymin':33.5,'ymax':35.0,'dx':0.5,'dy':0.5,'ny':4,'nx':5}
    data = np.arange(0,20).reshape((4,5))
    grid = Grid2D(data,geodict)
    
    yestestgeovalues = sampleShapes(shapes,YesTestPoints,'value')
    yestestgridvalues = sampleFromGrid(grid,YesTestPoints)
    
    notestgeovalues = sampleShapes(shapes,NoTestPoints,'value')
    notestgridvalues = sampleFromGrid(grid,YesTestPoints)

    TestFrame = pd.DataFrame()
    yescov = np.ones_like(yestestgeovalues[:,0])
    nocov = np.ones_like(notestgeovalues[:,0])
    yeslat = yestestgeovalues[:,1]
    yeslon = yestestgeovalues[:,0]
    nolat = notestgeovalues[:,1]
    nolon = notestgeovalues[:,0]
    TestFrame['Latitude'] = np.vstack((yeslat,nolat))
    TestFrame['Longitude'] = np.vstack((yeslat,nolat))
    TestFrame['Coverage'] = np.vstack((yescov,nocov))
    TestFrame['Geology'] = np.vstack((yestestgeovalues,notestgeovalues))
    
if __name__ == '__main__':
    _test_sample()
    
            
    
   
    
        
                
        
        
