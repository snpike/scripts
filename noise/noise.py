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

trigX = [(3*j) + 1  for j in range(11)]
trigY = [(3*j) + 1  for j in range(11)]

channel = []
channelMap = [[[] for i in range(32)] for j in range(32)]
rawx = []
rawy = []
for i in np.arange(START, END):
	if (not np.isnan(data['PH'][i])) and (0 < data['PH'][i] < 5000):
		channel.append(data['PH'][i])
		channelMap[data['RAWX'][i]][data['RAWY'][i]].append(data['PH'][i])
		rawx.append(data['RAWX'][i])
		rawy.append(data['RAWY'][i])

'''
spectrum = np.histogram(channel, bins=int(np.ceil(np.max(channel))))
plt.figure()
plt.plot(range(len(spectrum[0])), spectrum[0])
plt.show()
plt.close()'''

countMap = [[len(channelMap[i][j]) for i in range(32)] for j in range(32)]
plt.figure()
plt.imshow(countMap)
plt.colorbar()
plt.show()
plt.close()

'''
trigSum = 0

FWHM = []
for x in trigX:
	for y in trigY:
		if(channelMap[x][y]):
			trigSum += 1
			#print(x)
			#print(y)
			tempSpec = np.histogram(channelMap[x][y], bins=int(np.ceil(np.max(channelMap[x][y]))))
			centroid = np.argmax(tempSpec[0])
			fit_channels = np.arange(centroid-100, centroid + 250)
			g_init = models.Gaussian1D(amplitude=tempSpec[0][centroid], mean=centroid, stddev = 75)
			fit_g = fitting.LevMarLSQFitter()
			g = fit_g(g_init, range(len(tempSpec[0])), tempSpec[0])
			FWHM.append(2*np.sqrt(2*np.log(2))*g.stddev)
			#plt.plot(range(len(tempSpec[0])), tempSpec[0])
			#plt.plot(range(len(tempSpec[0])), g(range(len(tempSpec[0]))))
			#plt.show()
			#plt.close()
			
print(trigSum)	

FWHM_hist = np.histogram(FWHM, bins = 20)
plt.figure()
plt.plot(range(len(FWHM_hist[0])), FWHM_hist[0])
plt.show()
plt.close()
'''

noiseHist = np.histogram(np.array(countMap).flatten(), bins = np.array(1,15))
plt.figure()
plt.step(noiseHist[1][:-1], noiseHist[0], where='post')
plt.ylabel('Pixels')
plt.xlabel('Counts')
plt.show()
plt.close()

bins = np.arange(1,5000)
spectrum = np.histogram(data['PH'][START:END], bins = bins)
plt.plot(spectrum[1][:-1], spectrum[0])
plt.xlabel('Channel')
plt.ylabel('Counts')
plt.show()
plt.close()