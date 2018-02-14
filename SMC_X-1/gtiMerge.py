from astropy.io import fits
import numpy as np
import datetime

events = fits.open('/Users/sean/Desktop/SMC_X-1_filt/pipe_out003/nu10002013003A01_cl.evt', mode='update')
usrgtifile = fits.open('/Users/sean/Desktop/SMC_X-1_filt/pipe_out003/xselect_filter_15.0-40.0_fpma_gti.fits')

tref = usrgtifile[0].header['TIMEZERO']
usrgti = usrgtifile[1].data
events_gti = events[2].data

gti_filt = []
for i in range(len(events_gti)):
	start, stop = events_gti[i]
	for x in usrgti:
		x0, x1 = x
		x0 += tref
		x1 += tref
		if x0<=start<x1 and x0<stop<=x1:
			gti_filt.append([start, stop])
		if x0<=start<x1 and (not x0<stop<=x1):
			gti_filt.append([start, x1])
		if x0<stop<=x1 and (not x0<=start<x1):
			gti_filt.append([x0, stop])

rows = np.array(gti_filt)
columns = rows.T
col0 = fits.Column(name = 'START', format = 'f8', array = columns[0])
col1 = fits.Column(name = 'STOP', format = 'f8', array = columns[1])

cols = [col0, col1]

record = fits.FITS_rec.from_columns(cols)
events[2].data = record
events[0].header['history'] ='Filtered gti using xselect intensity filter (15.0 - 40.0 cts/s) xselect_filter_15.0-40.0_fpma_gti.fits on ' + str(datetime.datetime.now())
events.flush()