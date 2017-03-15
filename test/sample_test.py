#!/usr/bin/env python

#stdlib imports
import os.path
import sys
from collections import OrderedDict
import warnings

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

def sample_test():
    

if __name__ == '__main__':
    sample_test()
