from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from astropy.modeling import models, fitting
import os.path
import pickle

Am_line = 59.54

filepath = input('Please enter the filepath to the gamma flood data: ').strip()

slash = 0
i = 0
for char in filepath:
    if char == '/':
        slash = i
    i += 1

filename = filepath[slash + 1:]
detector = input('Please enter the detector ID: ').strip()
source = input('Please enter the source: ').strip()
etc = input('Please enter any other important information (temperature, voltage, etc.): ')

gainBool = os.path.exists('/disk/lif2/spike/detectorData/' + detector + '/' + filename[:-4] + 'quickgain.txt')

gain = np.zeros((34, 34))
gain[1:33][1:33] = pickle.load(open('/disk/lif2/spike/detectorData/' + detector + '/' + filename[:-4] + 'quickgain.txt', 'rb'))

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

countMap = [[np.sum(np.multiply(TOTmask, np.multiply(np.array(data['RAWX'][START:END])==i, np.array(data['RAWY'][START:END])==j))) for i in range(32)] for j in range(32)]

'''
for i in np.arange(START, END):
    if (not np.isnan(data['PH'][i])) and (0 < data['PH'][i] < maxchannel) and not data['STIM'][i]:
        countMap[data['RAWX'][i]][data['RAWY'][i]] += 1'''

plt.figure()
masked = np.ma.masked_values(countMap, 0.0)
current_cmap = mpl.cm.get_cmap()
current_cmap.set_bad(color='gray')
plt.imshow(masked)
c = plt.colorbar()
c.set_label('Counts')
#plt.title(detector + ' ' + source + ' Pixel Map ' + '(' + etc + ')')
plt.tight_layout()
plt.savefig('/disk/lif2/spike/detectorData/' + detector + '/figures/' + filename[:-4] + 'floodmap.eps')
#plt.show()
plt.close()

plt.figure()
plt.hist(np.array(countMap), bins = 100, range = (0, np.max(countMap) + 1), histtype = 'step')
plt.ylabel('Pixels')
plt.xlabel('Counts')
#plt.xticks(noiseHist[1][1:-1])
#plt.title(detector + ' ' + source + ' Count Histogram ' + '(' + etc + ')')
plt.tight_layout()
plt.savefig('/disk/lif2/spike/detectorData/' + detector + '/figures/' + filename[:-4] + 'gammahist.eps')
#plt.show()
plt.close()

bins = np.arange(1,maxchannel)
spectrum = np.histogram(data['PH'][START:END], bins = bins, range= (0, maxchannel))

centroid = np.argmax(spectrum[0][1000:]) + 1000
fit_channels = np.arange(centroid-100, centroid + 200)
g_init = models.Gaussian1D(amplitude=spectrum[0][centroid], mean=centroid, stddev = 75)
fit_g = fitting.LevMarLSQFitter()
g = fit_g(g_init, fit_channels, spectrum[0][fit_channels])
print(np.diag(fit_g.fit_info['param_cov']))
sigma_err = np.diag(fit_g.fit_info['param_cov'])[2]
fwhm_err = 2*np.sqrt(2*np.log(2))*sigma_err
mean_err = np.diag(fit_g.fit_info['param_cov'])[1]
frac_err = np.sqrt(np.square(fwhm_err) + np.square(g.fwhm*mean_err/g.mean))/g.mean
print(g.fwhm/g.mean)
print(frac_err)
print(Am_line * g.fwhm/g.mean)
print(frac_err * Am_line)
#plt.text(maxchannel*3/5, spectrum[0][centroid]*3/5, r'$\mathrm{FWHM}=$' + str(int(g.fwhm)) + r'$\pm$' + str(int(2*np.sqrt(2*np.log(2))*sigma_err)), fontsize=16)
#plt.text(maxchannel*3/5, spectrum[0][centroid]*3/5, r'$\mathrm{\frac{FWHM}{\mu}}=$' + str(int(round(100*g.fwhm/g.mean, 0))) + '%', fontsize=14)
plt.text(maxchannel*3/5, spectrum[0][centroid]*3/5, r'$\mathrm{FWHM}=$' + str(int(round(Am_line * 1000 * g.fwhm/g.mean, 0))) + r'$\pm$' + str(int(round(frac_err * Am_line*1000))) + ' eV', fontsize=13)

plt.plot(spectrum[1][:-1], spectrum[0], label = r'${}^{241}{\rm Am}$')
plt.plot(fit_channels, g(fit_channels), label = 'Gaussian fit')
plt.xlabel('Channel')
plt.ylabel('Counts')
plt.legend()

#plt.title(detector + ' ' + source + ' Spectrum ' + '(' + etc + ')')
plt.tight_layout()
plt.savefig('/disk/lif2/spike/detectorData/' + detector + '/figures/' + filename[:-4] + 'gammaspec.eps')
#plt.show()
plt.close()