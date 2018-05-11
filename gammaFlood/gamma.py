from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from astropy.modeling import models, fitting

filepath = input('Please enter the filepath to the gamma flood data: ')

slash = 0
i = 0
for char in filepath:
    if char == '/':
        slash = i
    i += 1

filename = filepath[slash + 1:]
detector = input('Please enter the detector ID: ')
source = input('Please enter the source: ')
etc = input('Please enter any other important information (temperature, voltage, etc.): ')

file = fits.open(filepath)

data = file[1].data

mask = data['TEMP'] > -20

START = np.argmax(mask)
END = len(mask) - np.argmax(mask[::-1])

trigX = [(3*j) + 1  for j in range(11)]
trigY = [(3*j) + 1  for j in range(11)]

maxchannel = 10000

PHmask = np.multiply(0 < np.array(data['PH'][START:END]),np.array(data['PH'][START:END]) < maxchannel)
STIMmask = np.array(data['STIM'][START:END])==0
TOTmask = np.multiply(PHmask, STIMmask)

countMap = [[np.sum(np.multiply(TOTmask, np.multiply(np.array(data['RAWX'])==i, np.array(data['RAWY'])==j))) for i in range(32)] for j in range(32)]

'''
for i in np.arange(START, END):
    if (not np.isnan(data['PH'][i])) and (0 < data['PH'][i] < maxchannel) and not data['STIM'][i]:
        countMap[data['RAWX'][i]][data['RAWY'][i]] += 1'''

plt.figure()
plt.imshow(countMap)
c = plt.colorbar()
c.set_label('Counts')
plt.title(detector + ' ' + source + ' Pixel Map ' + '(' + etc + ')')
plt.tight_layout()
plt.savefig('/disk/lif2/spike/detectorData/' + detector + '/figures/' + filename[:-4] + 'floodmap.eps')
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

noiseHist = np.histogram(np.array(countMap).flatten(), bins = np.arange(0,np.max(np.array(countMap).flatten()) + 3))
plt.figure()
plt.step(noiseHist[1][:-1], noiseHist[0], where='mid')
plt.xlim(0.5, noiseHist[1][-1])
plt.ylim(-1,np.max(noiseHist[0][1:]) + 1)
plt.ylabel('Pixels')
plt.xlabel('Counts')
#plt.xticks(noiseHist[1][1:-1])
plt.title(detector + ' ' + source + ' Count Histogram ' + '(' + etc + ')')
plt.tight_layout()
plt.savefig('/disk/lif2/spike/detectorData/' + detector + '/figures/' + filename[:-4] + 'gammahist.eps')
plt.show()
plt.close()

bins = np.arange(1,maxchannel)
spectrum = np.histogram(data['PH'][START:END], bins = bins, range= (0, maxchannel))

centroid = np.argmax(spectrum[0][1000:]) + 1000
fit_channels = np.arange(centroid-200, centroid + 250)
g_init = models.Gaussian1D(amplitude=spectrum[0][centroid], mean=centroid, stddev = 75)
fit_g = fitting.LevMarLSQFitter()
g = fit_g(g_init, range(len(spectrum[0])), spectrum[0])
plt.text(maxchannel*3/4, spectrum[0][centroid]*4/5, r'$\mathrm{FWHM}=$' + str(int(2*np.sqrt(2*np.log(2))*g.stddev)), fontsize=16)

plt.plot(spectrum[1][:-1], spectrum[0])
plt.plot(fit_channels, g(fit_channels))
plt.xlabel('Channel')
plt.ylabel('Counts')
plt.title(detector + ' ' + source + ' Spectrum ' + '(' + etc + ')')
plt.tight_layout()
plt.savefig('/disk/lif2/spike/detectorData/' + detector + '/figures/' + filename[:-4] + 'gammaspec.eps')
plt.show()
plt.close()