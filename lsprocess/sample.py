#!/usr/bin/env python

#stdlib imports
import os.path
from functools import partial
from collections import OrderedDict

import numpy as np
import fiona
from shapely.geometry import Polygon,shape,MultiPoint,Point
import matplotlib.pyplot as plt
import pyproj
import pandas as pd
from shapely.ops import transform
from grid.grid2d import Grid2D
from grid.gmt import GMTGrid
from grid.gdal import GDALGrid
from grid.shake import ShakeGrid

def getFileType(filename):
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
    latmiddle = ymin + (ymax-ymin)/2.0
    lonmiddle = xmin + (xmax-xmin)/2.0
    projstr = '+proj=ortho +lat_0=%.4f +lon_0=%.4f +x_0=0.0 +y_0=0.0' % (latmiddle,lonmiddle)
    proj = pyproj.Proj(projparams=projstr)
    project = partial(
        pyproj.transform,
        pyproj.Proj(projparams='+proj=latlong'),
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
            yespoints.append(pshape['geometry']['coordinates'])
            
    return (np.array(yespoints),nrows,ncols,xvar,yvar,idx)

def sampleFromFile(shapefile,dx=10.0,nmax=None,testPercent=1.0,classBalance=None,extent=None,Nsamp=100):
    """
    shapefile - path to shapefile, presumed to be geographic
    dx - resolution of sampling in X and Y
    nmax - if not None, maximum allowed number of mesh points in X and Y together (nrows*ncols).  Overrides dx.
    testPercent - says how many samples to put into output data frame.
    classBalance - If None, uses CB of data (raise exception for points).  If specified, is the fraction of "yes" pixels, where "yes"+"no" = 1
    extent - if none, use the bounding box of the data in the shapefile. 
    """
    #read the shapes in from the file
    f = fiona.collection(shapefile,'r')
    shapes = list(f)
    bounds = f.bounds
    
    f.close()

    return sampleFromShapes(shapes,bounds,dx=dx,nmax=nmax,testPercent=testPercent,classBalance=classBalance,extent=extent,Nsamp=Nsamp)

def sampleYes(array,N):
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

def sampleFromShapes(shapes,bounds,dx=10.0,nmax=None,testPercent=1.0,classBalance=None,extent=None,Nsamp=100):
    xmin,ymin,xmax,ymax = bounds
    shptype = shapes[0]['geometry']['type']
    if shptype not in ['Point','Polygon']:
        raise Exception('Only polygon and point data types supported!')
    
    if shptype != 'Polygon':
        if classBalance is None:
            raise Exception('Class balance *must* be selected when input data are points.')

    #Get the shapes projected into an orthographic projection centered on the data 
    pshapes,proj = getProjectedShapes(shapes,xmin,xmax,ymin,ymax)

    #get the "yes" sample points
    yespoints,nrows,ncols,xvar,yvar,yesidx = getYesPoints(pshapes,proj,dx,nmax)

    #what is the class balance (assuming polygons)
    if classBalance is None:
        classBalance = getClassBalance(pshapes,bounds,proj)

    if extent is None: #we're using the default bounding box of the coverage data
        Nmesh = nrows*ncols
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

def projectBack(points,proj):
    mpoints = MultiPoint(points)
    project = partial(
        pyproj.transform,
        proj,
        pyproj.Proj(projparams='+proj=latlong'))
    gmpoints = transform(project, mpoints)
    coords = []
    for point in gmpoints.geoms:
        x,y = point.coords[0]
        coords.append((x,y))
    coords = np.array(coords)
    return coords

def plotPoints(shapes,YesTestPoints,YesTrainPoints,NoTestPoints,NoTrainPoints,filename):
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
    xmin,ymin,xmax,ymax = bounds
    bpoly = Polygon([(xmin,ymax),
                     (xmax,ymax),
                     (xmax,ymin),
                     (xmin,ymin)])
    project = partial(
        pyproj.transform,
        pyproj.Proj(projparams='+proj=latlong'),
        proj)
    bpolyproj = transform(project,bpoly)
    totalarea = bpolyproj.area
    polyarea = 0
    for pshape in pshapes:
        polyarea += pshape.area

    return polyarea/totalarea

def sampleShapeFile(shapefile,xypoints,attribute):
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
    xmin = xmin - fdict['xdim']*3
    xmax = xmax + fdict['xdim']*3
    ymin = ymin - fdict['ydim']*3
    ymax = ymax + fdict['ydim']*3
    sdict = {'xmin':xmin,'xmax':xmax,'ymin':ymin,'ymax':ymax}
    if gridtype == 'gmt':
        grid = GMTGrid.load(gridfile,samplegeodict=sdict,resample=False,method=method,doPadding=True)
    else:
        grid = GDALGrid.load(gridfile,samplegeodict=sdict,resample=False,method=method,doPadding=True)

    return sampleFromGrid(grid,xypoints)

def sampleFromGrid(grid,xypoints):
    samples = []
    for lon,lat in xypoints:
        sample = grid.getValue(lat,lon)
        samples.append(sample)

    return np.array(samples)

def sampleFromShakeMap(shakefile,layer,xypoints):
    shakegrid = ShakeGrid.load(shakefile)
    return sampleFromMultiGrid(shakegrid,layer,points)

def sampleFromMultiGrid(multigrid,layer,xypoints):
    if layer not in multigrid.getLayerNames():
        raise Exception('Layer %s not found in grid' % layer)
    hazgrid = multigrid.getLayer(layer)
    return sampleFromGrid(hazgrid,xypoints)

def getDataFrames(sampleparams,shakeparams,predictors,outparams):
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

    yestest,yestrain,notest,notrain,xvar,yvar,pshapes,proj = sampleFromFile(coverage,dx=dx,nmax=nmax,testPercent=testpercent,classBalance=cb,extent=extent,Nsamp=nsamp)

    traincolumns = OrderedDict()
    testcolumns = OrderedDict()
    traincolumns['lat'] = np.concatenate((yestrain[:,1],notrain[:,1]))
    traincolumns['lon'] = np.concatenate((yestrain[:,0],notrain[:,0]))
    traincolumns['coverage'] = np.concatenate((np.ones_like(yestrain[:,1]),np.ones_like(notrain[:,1])))
    testcolumns['lat'] = np.concatenate((yestest[:,1],notest[:,1]))
    testcolumns['lon'] = np.concatenate((yestest[:,0],notest[:,0]))
    testcolumns['coverage'] = np.concatenate((np.ones_like(yestest[:,1]),np.ones_like(notest[:,1])))
    
    
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
            yes_test_samples = sampleGridFile(predfile,yestest,method=method)
            no_test_samples = sampleGridFile(predfile,notest,method=method)
            yes_train_samples = sampleGridFile(predfile,yestrain,method=method)
            no_train_samples = sampleGridFile(predfile,notrain,method=method)
            testcolumns[predname] = np.concatenate((yes_test_samples,no_test_samples))
            traincolumns[predname] = np.concatenate((yes_train_samples,no_train_samples))
        else:
            continue #attribute or sampling method key

    #sample the shakemap
    layers = ['mmi','pga','pgv','psa03','psa10','psa30']
    shakegrid = ShakeGrid.load(shakeparams['shakemap'])
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
    geodict = {'xmin':33.5,'xmax':35.5,'ymin':33.5,'ymax':35.0,'xdim':0.5,'ydim':0.5,'nrows':4,'ncols':5}
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
    
            
    
   
    
        
                
        
        
