from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from astropy.modeling import models, fitting

filepath = input('Please enter the filepath to the noise data:')

file = fits.open(filepath)

data = file[1].data

mask = data['TEMP'] > -20

START = np.argmax(mask)
END = len(mask) - np.argmax(mask[::-1])

trigX = [3*j for j in range(5)] + [3*j + 16 for j in range(6)]
trigY = [3*j for j in range(6)] + [3*j + 19 for j in range(5)]

channel = []
channelMap = [[[] for i in range(32)] for j in range(32)]
rawx = []
rawy = []
for i in np.arange(START, END):
	if (not np.isnan(data['PH'][i])) and (0 < data['PH'][i] < 50000):
		channel.append(data['PH'][i])
		channelMap[data['RAWX'][i]][data['RAWY'][i]].append(data['PH'][i])
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

trigSum = 0

for x in trigX:
	for y in trigY:
		if(channelMap[x][y]):
			trigSum += 1
			#print(x)
			#print(y)
			tempSpec = np.histogram(channelMap[x][y], bins=int(np.ceil(np.max(channelMap[x][y]))))
			centroid = np.argmax(tempSpec[0])
			#fit_channels = np.arange(centroid-100, centroid + 250)
			g_init = models.Gaussian1D(amplitude=tempSpec[0][centroid], mean=centroid, stddev = 75)
			fit_g = fitting.LevMarLSQFitter()
			g = fit_g(g_init, range(len(tempSpec[0])), tempSpec[0])
			#plt.plot(range(len(tempSpec[0])), tempSpec[0])
			#plt.plot(range(len(tempSpec[0])), g(range(len(tempSpec[0]))))
			#plt.show()
			#plt.close()
			
print(trigSum)	

'''
noiseHist = np.histogram(pixelmap[0].flatten(), bins = 20)
plt.figure()
plt.step(noiseHist[1][:-1], noiseHist[0], where='post')
plt.show()
plt.close()'''