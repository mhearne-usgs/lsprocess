from distutils.core import setup

setup(name='lsprocess',
      version='0.1dev',
      description='NEIC Landslide Data Processing Tools',
      author='Mike Hearne',
      author_email='mhearne@usgs.gov',
      url='',
      install_requires=['numpy'],
      scripts = ['lsprocess.py','xy2shp.py'],
)
