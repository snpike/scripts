#! /usr/bin/env python
# Given two positions, give offsets and PA.

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
import argparse

parser = argparse.ArgumentParser(description='Calculate offsets from position 2 to position 1, and the PA of a slit covering both. \
	-d means RA, Dec is in degrees, otherwise sex is assumed. If you want to input negative declinations, preceed with -- (otherwise python interprets this as an argument).')
parser.add_argument('position1', nargs=2, type=str, help='RA, DEC of target; ddd.ddddd or hh:mm:ss.sss')
parser.add_argument('position2', nargs=2, type=str, help='RA, DEC of offset star or second target to be placed in slit')
parser.add_argument('-d', '--deg', action='store_true', default='False', help='Input in decimal degrees')

args = parser.parse_args()

if args.deg == True:
    pos1 = SkyCoord(float(args.position1[0]),float(args.position1[1]), unit='deg')
    pos2 = SkyCoord(float(args.position2[0]),float(args.position2[1]), unit='deg')
else:
    pos1 = SkyCoord(args.position1[0],args.position1[1], unit=(u.hourangle, u.deg)) 
    pos2 = SkyCoord(args.position2[0],args.position2[1], unit=(u.hourangle, u.deg)) 

pa = pos1.position_angle(pos2).degree
offra = pos2.spherical_offsets_to(pos1)[0].arcsec
offdec = pos2.spherical_offsets_to(pos1)[1].arcsec

print( "Position angle:", pa, "degrees")
print( "Offset from offset star to target:")
if offra > 0:
    print(offra, "arcsec East")
else:
    print( abs(offra), "arcsec West")
if offdec > 0:
    print( offdec, "arcsec North")
else:
    print( abs(offdec), "arcsec South")