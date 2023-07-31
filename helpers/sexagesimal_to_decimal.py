#!/usr/bin/env python
# coding: utf-8

import numpy as np

ra = input("RA: ").split()
while len(ra) !=3:
	print("Invalid RA input. Try again.")
	ra = input("RA: ").split()
dec = input('Dec: ').split()
while len(dec) !=3:
	print("Invalid RA input. Try again.")
	dec = input("RA: ").split()


ra_decimal = float(ra[0])*(360./24.) + float(ra[1])*(360./(24.*60)) + float(ra[2])*(360./(24.*3600.))

dec_decimal = float(dec[0]) + float(dec[1])/60. + float(dec[2])/3600.

print("RA in decimal: ", str(ra_decimal))
print("Dec in decimal: ", str(dec_decimal))