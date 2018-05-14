from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from astropy.modeling import models, fitting

filepath = input('Please enter the filepath to the noise data: ')

slash = 0
i = 0
for char in filepath:
	if char == '/':
		slash = i
	i += 1

filename = filepath[slash + 1:]
detector = input('Please enter the detector ID: ')
etc = input('Please enter any other important information (temperature, voltage, etc.): ')

test = 'noise'

file = fits.open(filepath)

data = file[1].data

mask = data['TEMP'] > -20

START = np.argmax(mask)
END = len(mask) - np.argmax(mask[::-1])

trigX = [(3*j) + 1  for j in range(11)]
trigY = [(3*j) + 1  for j in range(11)]


maxchannel = 2000
bins = np.arange(1,maxchannel)

channel = []
channelMap = [[[] for i in range(32)] for j in range(32)]
rawx = []
rawy = []
for i in np.arange(START, END):
	if data['UP'][i]:
		channelMap[data['RAWX'][i]][data['RAWY'][i]].append(data['PH'][i])
		rawx.append(data['RAWX'][i])
		rawy.append(data['RAWY'][i])

countMap = [[len(channelMap[i][j]) for i in range(32)] for j in range(32)]
plt.figure()
plt.imshow(countMap)
c = plt.colorbar()
c.set_label('Counts')
plt.title(detector + ' ' + test + ' Pixel Map ' + '(' + etc + ')')
plt.tight_layout()
plt.savefig('/disk/lif2/spike/detectorData/' + detector + '/figures/' + filename[:-4] + 'pixmap.eps')
#plt.show()
plt.close()

FWHM = []
FWHM_map = [[np.nan for i in range(32)] for j in range(32)]
for x in trigX:
	for y in trigY:
		if(channelMap[x][y]):
			tempSpec = np.histogram(channelMap[x][y], bins=bins, range = (0,maxchannel*2))
			centroid = np.argmax(tempSpec[0])
			fit_channels = np.arange(centroid-100, centroid + 250)
			g_init = models.Gaussian1D(amplitude=tempSpec[0][centroid], mean=centroid, stddev = 75)
			fit_g = fitting.LevMarLSQFitter()
			g = fit_g(g_init, range(len(tempSpec[0])), tempSpec[0])
			FWHM.append(2*np.sqrt(2*np.log(2))*g.stddev)
			if 2*np.sqrt(2*np.log(2))*g.stddev < 1000:
				FWHM_map[x][y] = 2*np.sqrt(2*np.log(2))*g.stddev
			plt.step(tempSpec[1][:-1], tempSpec[0], where='mid')
			plt.plot(fit_channels, g(fit_channels))
			plt.ylabel('Counts')
			plt.xlabel('Channel')
			plt.tight_layout()
			plt.savefig('/disk/lif2/spike/detectorData/' + detector + '/figures/pixel_figs/' + filename[:-4] + 'x' + str(x) + 'y' + str(y) + '_spec.eps')
			plt.close()

FWHM_hist = np.histogram(FWHM, bins = 20, range = (0, 500))
plt.figure()
plt.step(FWHM_hist[1][:-1], FWHM_hist[0], where='mid')
plt.ylabel('Pixels')
plt.xlabel('FWHM (channels)')
plt.title(detector + ' ' + test + ' FWHM Histogram ' + '(' + etc + ')')
plt.tight_layout()
plt.savefig('/disk/lif2/spike/detectorData/' + detector + '/figures/' + filename[:-4] + 'FWHMhist.eps')
#plt.show()
plt.close()

plt.figure()
current_cmap = matplotlib.cm.get_cmap()
current_cmap.set_bad(color='gray')
plt.imshow(FWHM_map)
c = plt.colorbar()
c.set_label('FWHM (channels)')
plt.title(detector + ' ' + test + ' FWHM Map ' + '(' + etc + ')')
plt.tight_layout()
plt.savefig('/disk/lif2/spike/detectorData/' + detector + '/figures/' + filename[:-4] + 'FWHMmap.eps')
plt.show()
plt.close()


noiseHist = np.histogram(np.array(countMap).flatten(), bins = np.arange(0,np.max(np.array(countMap).flatten()) + 3))
plt.figure()
plt.step(noiseHist[1][:-1], noiseHist[0], where='mid')
plt.xlim(0.5, noiseHist[1][-1])
plt.ylim(-1,np.max(noiseHist[0][1:]) + 1)
plt.ylabel('Pixels')
plt.xlabel('Counts')
#plt.xticks(noiseHist[1][1:-1])
plt.title(detector + ' ' + test + ' Count Histogram ' + '(' + etc + ')')
plt.tight_layout()
plt.savefig('/disk/lif2/spike/detectorData/' + detector + '/figures/' + filename[:-4] + 'pixhist.eps')
#plt.show()
plt.close()

spectrum = np.histogram(data['PH'][START:END], bins = bins, range= (0, maxchannel*2))
plt.plot(spectrum[1][:-1], spectrum[0])
plt.xlabel('Channel')
plt.ylabel('Counts')
plt.title(detector + ' ' + test + ' Spectrum ' + '(' + etc + ')')
plt.tight_layout()
plt.savefig('/disk/lif2/spike/detectorData/' + detector + '/figures/' + filename[:-4] + 'spec.eps')
#plt.show()
plt.close()