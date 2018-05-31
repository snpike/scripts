from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from astropy.modeling import models, fitting
import pickle

# From http://www.nndc.bnl.gov/nudat2/indx_dec.jsp
lines = {'Am241': [59.54, r'${}^{241}{\rm Am}$'], 'Co57': [122.06, r'${}^{57}{\rm Co}$']}

filepath = input('Please enter the filepath to the gamma flood data: ').strip()

slash = 0
i = 0
for char in filepath:
    if char == '/':
        slash = i
    i += 1

filename = filepath[slash + 1:]
detector = input('Please enter the detector ID: ').strip()
source = input('Please enter the source (Am241 or Co57): ').strip()
while source not in lines:
	source = input('Please enter the source (Am241 or Co57): ').strip()
etc = input('Please enter any other important information (temperature, voltage, etc.): ')

file = fits.open(filepath)

data = file[1].data

mask = data['TEMP'] > -20

START = np.argmax(mask)
END = len(mask) - np.argmax(mask[::-1])

maxchannel = 10000
bins = np.arange(1,maxchannel)

gain = np.zeros((32,32))

for x in range(32):
	for y in range(32):

		plt.figure()

		channel = data.field('PH')[START:END][np.nonzero(np.multiply((data.field('RAWX')[START:END] == x), (data.field('RAWY')[START:END] == y)))]

		if len(channel):
			spectrum = np.histogram(channel, bins=bins, range = (0, maxchannel))
			centroid = np.argmax(spectrum[0][3000:]) + 3000
			fit_channels = np.arange(centroid-100, centroid + 200)
			g_init = models.Gaussian1D(amplitude=spectrum[0][centroid], mean=centroid, stddev = 75)
			fit_g = fitting.LevMarLSQFitter()
			g = fit_g(g_init, fit_channels, spectrum[0][fit_channels])
			if fit_g.fit_info['param_cov'] is not None:
				sigma_err = np.diag(fit_g.fit_info['param_cov'])[2]
				fwhm_err = 2*np.sqrt(2*np.log(2))*sigma_err
				mean_err = np.diag(fit_g.fit_info['param_cov'])[1]
				frac_err = np.sqrt(np.square(fwhm_err) + np.square(g.fwhm*mean_err/g.mean))/g.mean
				plt.text(maxchannel*3/5, spectrum[0][centroid]*3/5, r'$\mathrm{FWHM}=$' + str(int(round(lines[source][0] * 1000 * g.fwhm/g.mean, 0))) + r'$\pm$' + str(int(round(frac_err * lines[source][0]*1000))) + ' eV', fontsize=13)
				gain[y][x] = lines[source][0]/g.mean

				plt.hist(np.multiply(channel,lines[source][0]/g.mean), bins=np.multiply(bins,lines[source][0]/g.mean), range = (0, maxchannel * lines[source][0]/g.mean), histtype='step')
				plt.plot(np.multiply(fit_channels,lines[source][0]/g.mean), g(fit_channels), label = 'Gaussian fit')
				plt.ylabel('Counts')
				plt.xlabel('Energy')
				plt.legend()

		plt.tight_layout()		
		plt.savefig('/disk/lif2/spike/detectorData/' + detector + '/figures/pixel_figs/' + filename[:-4] + '_x' + str(x) + '_y' + str(y) + '_gammaspec.eps')
		plt.close()

# interpolate gain for pixels where fit was unsuccessful
interpgain = []
newgain = np.zeros((34,34))
newgain[1:33, 1:33] = gain
empty = np.nonzero(gain == 0)
for i in range(len(empty)):
	temp = newgain[empty[i][0]:empty[i][0]+3, empty[i][1]:empty[i][1]+3]
	interpgain.append(np.sum(temp)/np.count_nonzero(temp))

for i in range(len(empty)):
	gain[empty[i]] = interpgain[i]

pickle.dump(gain, open('/disk/lif2/spike/detectorData/' + detector + '/' + filename[:-4] + 'quickgain.txt', 'wb'))

plt.figure()
plt.imshow(gain)
c = plt.colorbar()
c.set_label('keV/Channels')
#plt.title(detector + ' ' + source + ' Pixel Map ' + '(' + etc + ')')
plt.tight_layout()
#plt.show()
plt.savefig('/disk/lif2/spike/detectorData/' + detector + '/figures/' + filename[:-4] + 'gainmap.eps')
plt.close()