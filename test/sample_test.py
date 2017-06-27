#!/usr/bin/env python

#stdlib imports
import os.path
import sys
from collections import OrderedDict
import warnings
import glob

#hack the path so that I can debug these functions if I need to
homedir = os.path.dirname(os.path.abspath(__file__)) #where is this script?
libdir = os.path.abspath(os.path.join(homedir,'..'))
sys.path.insert(0,libdir) #put this at the front of the system path, ignoring any installed sampling stuff

import numpy as np
import fiona
from shapely.geometry import Polygon,shape,MultiPoint,Point
import os.path
import matplotlib.pyplot as plt
from functools import partial
import pyproj
from shapely.ops import transform

from lsprocess import sample

def sample_test_points():
    datadir = os.path.join(homedir,'data')
    coverage = os.path.join(datadir,'Keefer_and_Manson_1998.geojson')
    shakemap = os.path.join(datadir,'lomaprieta_grid.xml')
    popfile = os.path.join(datadir,'lomaprieta_pop.bil')
    params = {}
    params['coverage'] = coverage
    params['dx'] = 50
    params['cb'] = 0.5
    params['nmax'] = 150000000
    params['nsamp'] = 1000
    params['touch_center'] = True
    params['testpercent'] = 0
    params['h1'] = 85
    params['h2'] = 395
    params['extent'] = (-125.380000,-119.380000,34.646000,39.434000)
    shakeparams = {'shakemap':shakemap}
    predictors = {'population':popfile,'population_sampling':'nearest'}
    outparams = {'folder':os.path.expanduser('~'),'basename':'lomaprieta'}
    training,testing = sample.getDataFrames(params,shakeparams,predictors,outparams)

def sample_test_polygons():
    datadir = os.path.join(homedir,'data')
    coverage = os.path.join(datadir,'Harp_and_Jibson_1995.geojson')
    shakemap = os.path.join(datadir,'northridge.xml')
    popfile = os.path.join(datadir,'northridge_gpw.flt')
    params = {}
    params['coverage'] = coverage
    params['dx'] = 10
    params['cb'] = 0.5
    params['nmax'] = 10000
    params['nsamp'] = 100000
    params['touch_center'] = False
    params['testpercent'] = 0.5
    params['h1'] = 10
    params['h2'] = 100
    #params['extent'] = (-121.046000,-116.046000,32.143500,36.278500)
    params['extent'] = None
    shakeparams = {'shakemap':shakemap}
    predictors = {'population':popfile,'population_sampling':'nearest'}
    outparams = {'folder':os.path.expanduser('~'),'basename':'northridge'}
    training,testing = sample.getDataFrames(params,shakeparams,predictors,outparams)

if __name__ == '__main__':
    sample_test_points()
    sample_test_polygons()
    
