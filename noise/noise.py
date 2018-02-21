from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np

filepath = input('Please enter the filepath to the noise data:')

file = fits.open(filepath)

data = file[1].data

mask = data['TEMP'] > -20

START = np.argmax(mask)
END = len(mask) - np.argmax(mask[::-1])

channel = []
rawx = []
rawy = []
for i in np.arange(START, END):
	if (not np.isnan(data['PH'][i])) and (0 < data['PH'][i] < 50000):
		channel.append(data['PH'][i])
		rawx.append(data['RAWX'][i])
		rawy.append(data['RAWY'][i])


spectrum = np.histogram(channel, bins=int(np.ceil(np.max(channel))))
plt.figure()
plt.plot(range(len(spectrum[0])), spectrum[0])
plt.show()
plt.close()

pixelmap = np.histogram2d(rawx, rawy, bins = (32,32))
plt.figure()
plt.imshow(pixelmap[0])
plt.colorbar()
plt.show()
plt.close()