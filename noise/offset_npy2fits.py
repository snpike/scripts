from astropy.io import fits
import numpy as np
import os.path

filepath = input('Enter the path to the .npy offset file to convert:').strip()
npy_data = np.load(filepath)

while npy_data.shape != (32,32,16):
	print('The shape of the offset array should be (32, 32, 16). Instead got ' + str(offsets.shape) + '. Try again.')
	filepath = input('Enter the path to the .npy offset file to convert:').strip()
	npy_data = np.load(filepath)

slash = 0
i = 0
for char in filepath:
    if char == '/':
        slash = i
    i += 1

filename = filepath[slash + 1:]

detector = input('Please enter the detector ID: ').strip()

rawx = []
rawy = []
offset = []

for row in range(32):
	for col in range(32):
		rawx.append(col)
		rawy.append(row)
		offset.append(npy_data[row][col])

col1 = fits.Column(name = 'RAWX', format = '1B', array = rawx)
col2 = fits.Column(name = 'RAWY', format = '1B', array = rawy)
col3 = fits.Column(name = 'OFFSET', format = '16E', array = offset)

coldefs = fits.ColDefs([col1, col2, col3])

hdu = fits.BinTableHDU.from_columns(coldefs)
hdu1 = fits.BinTableHDU.from_columns(coldefs)
hdu2 = fits.BinTableHDU.from_columns(coldefs)
hdu3 = fits.BinTableHDU.from_columns(coldefs)

hdu.append(hdu1)
hdu.append(hdu2)
hdu.append(hdu3)

hdu.writeto(filepath[:-4] + '.fits')