#!/bin/bash

VENV=pysample
PYVER=3.5

DEPARRAY=(numpy scipy matplotlib rasterio pandas shapely h5py gdal pytest pytest pytest-cov pytest-mpl jupyter ipython fiona pyproj)

#turn off whatever other virtual environment user might be in
source deactivate
    
#remove any previous virtual environments called lsprocess
CWD=`pwd`
cd $HOME;
conda remove --name $VENV --all -y
cd $CWD
    
#create a new virtual environment called $VENV with the below list of dependencies installed into it
conda create --name $VENV --yes --channel conda-forge python=3.5 ${DEPARRAY[*]} -y

#activate the new environment
source activate $VENV

#install some items separately
#conda install -y sqlalchemy #at the time of this writing, this is v1.0, and I want v1.1
conda install -y psutil

#do pip installs of those things that are not available via conda.
#do pip installs of those things that are not available via conda.
curl --max-time 60 --retry 3 -L https://github.com/gem/oq-engine/archive/master.zip -o openquake.zip
pip -v install --no-deps openquake.zip
rm openquake.zip

curl --max-time 60 --retry 3 -L https://github.com/usgs/MapIO/archive/0.6.1.zip -o mapio.zip
pip install mapio.zip
rm mapio.zip

#tell the user they have to activate this environment
echo "Type 'source activate ${VENV}' to use this new virtual environment."
