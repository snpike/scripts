from astropy.io import fits

gti = input('Enter the GTI .txt file:').split('.')[0]
tstart = []
tstop = []
with open(gti + ".txt") as file:
    for line in file:
        splitline = line.split()
        tstart.append(float(splitline[0]))
        tstop.append(float(splitline[1]))

primary_hdr = fits.Header()
primary_hdr['HDUCLASS'] = 'OGIP'
primary_hdr['HDUDOC'] = 'OGIP/92-009'
primary_hdu = fits.PrimaryHDU(header=primary_hdr)
gti_hdr = fits.Header()
gti_hdr['HDUCLASS'] = 'OGIP'
gti_hdr['EXTNAME'] = 'GTI'
gti_hdr['TELESCOP'] = 'NuSTAR'

c1 = fits.Column(name='START', array=tstart, format='D')
c2 = fits.Column(name='STOP', array=tstop, format='D')
gti_hdu = fits.BinTableHDU.from_columns([c1, c2], header = gti_hdr)
hdul = fits.HDUList([primary_hdu, gti_hdu])
hdul.writeto(gti + '.fits', overwrite=True)
