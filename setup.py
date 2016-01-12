from distutils.core import setup

setup(name='lsprocess',
      version='0.2dev',
      description='NEIC Landslide/Liquefaction Data Processing Tools',
      author='Mike Hearne',
      author_email='mhearne@usgs.gov',
      url='',
      packages=['lsprocess'],
      scripts = ['secsample','subsample'],
)
