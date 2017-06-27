#!/usr/bin/env python

#stdlib imports
import argparse
import os.path

#third party imports
from mpl_toolkits.basemap import shapefile

def main(args):
    if not os.path.isfile(args.xyfile):
        raise Exception('%s is not a valid file!' % args.xyfile)
    fpath,fext = os.path.splitext(args.xyfile)
    #need to implement writer in shapefile object or just use joel lawheads code directly...
    writer = shapefile.Writer(shapefile.POINT)
    writer.autoBalance = 1
    writer.field('Id',fieldType='N')
    f = open(args.xyfile,'rt')
    if args.skipFirstLine:
        f.readline()
    idnum = 1
    for line in f.readlines():
        lonstr,latstr = line.split(',')
        lon = float(lonstr)
        lat = float(latstr)
        writer.point(lon,lat)
        writer.record(idnum)
        idnum += 1
    f.close()
    writer.save(fpath)

if __name__ == '__main__':
    desc = """Convert a GMT-style .xy file of points into a Shapefile.
    """
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("xyfile", type=str,
                        help="A GMT-style xy file")
    parser.add_argument("-s", "--skip-first-line", dest='skipFirstLine',
                        action="store_true",help="Skip first line")
    args = parser.parse_args()
    main(args)
    
