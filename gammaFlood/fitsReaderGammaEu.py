from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.io import readsav
import pickle

highELine = 86.545

def gauss(x, a, b, c):
	pdf =  a*np.exp(np.square(x-b)/((-2.0)*np.square(c)))
	#return pdf*cdf + k
	return pdf

FWHM = 2.0*np.sqrt(np.log(4.0))

### Apply the gain/offset corrections to the data then pickle it for later use

gainOffset =  readsav('/Users/sean/CdTe/CdTeData/h100_-10C_500V_gain.sav')['calstruct']
PHGainOffset = readsav('/Users/sean/CdTe/CdTeData/h100_-10C_500V_fit_2.sav')['calstruct3'].reshape(32,32)

tempgain = []
tempoffset = []

for i in range(len(gainOffset)):
	tempgain.append(gainOffset[i]['gain'])
	tempoffset.append(gainOffset[i]['offset'])

gain = np.zeros((34, 34))
offset = np.zeros((34, 34))

gain[1:33, 1:33] = np.array(tempgain).reshape(32, 32)
offset[1:33, 1:33] = np.array(tempoffset).reshape(32, 32)

del tempoffset
del tempgain
del gainOffset

file = fits.open('/Users/sean/CdTe/gammaFlood/20170310_H100_gamma_eu155_-10C.0V.fits')
data = file[1].data
file.close()
START = 0
i = 0
events = {'channel': [], 'rawx': [], 'rawy': [], 'grade': []}

# Skip the beginning
while data.field('TEMP')[i] < -50:
	i += 1
START = i

while i < len(data.field('PH')):
	rawx = data.field('RAWX')[i]
	rawy = data.field('RAWY')[i]
	grade = data.field('GRADE')[i]
	### Apply the gain and offset ONLY to pixels with positive counts
	temp = data.field('PH_COM')[i].reshape(3,3)
	if np.sum(temp) > 0:
		mask = (temp > 0).astype(int)
		channel = np.sum(np.add(np.multiply(np.multiply(mask, temp), gain[rawx:rawx + 3, rawy:rawy + 3]), np.multiply(mask, offset[rawx:rawx + 3, rawy:rawy + 3])))
		if grade < 5:
			channel = (channel * PHGainOffset[rawx][rawy]['gain'][grade]) + PHGainOffset[rawx][rawy]['offset'][grade]
		if (not np.isnan(channel)) and (channel < 3500):
			events['channel'].append(channel)
			events['rawx'].append(rawx)
			events['rawy'].append(rawy)
			events['grade'].append(grade)

	i += 1


pickle.dump(events, open('/Users/sean/CdTe/gammaFlood/20170310_H100_gamma_eu155_-10C_dump.0V.pkl', 'wb'))
'''
### Load the pickled data and do some plotting/fit the 59 keV line to a gaussian. Use this fit to convert from channels to energy.

events = pickle.load(open('/Users/sean/CdTe/gammaFlood/20170310_H100_gamma_eu155_-10C_dump.0V.pkl', 'rb'))

spectrum = np.histogram(events['channel'], bins = int(np.max(events['channel'])))

plt.plot(range(len(spectrum[0])), np.divide(spectrum[0], np.max(spectrum[0])), lw = 0.5)
plt.xlabel('Channel')
plt.ylabel('Normalized Counts')
plt.savefig('/Users/sean/CdTe/gammaFlood/images/20170310_H100_gamma_eu155_-10C_spectrum.0V.eps')
plt.close()


centroid = np.argmax(spectrum[0])
fitBounds = [centroid - 7, centroid + 15]
plotBounds = [fitBounds[0] - 50, fitBounds[1] + 42]
popt, pcov = curve_fit(gauss, range(fitBounds[0], fitBounds[1]), np.divide(spectrum[0], np.max(spectrum[0]))[fitBounds[0]:fitBounds[1]], p0 = [1.0, centroid, 1.0])
conversion = highELine/popt[1]

#stdev error
errors = np.sqrt(np.diag(pcov))

conversionErr = errors[1]*highELine/np.square(popt[1])

widthError = FWHM*np.sqrt(np.square(errors[2]*conversion) + np.square(conversionErr*popt[2]))


plt.plot(np.multiply(conversion,range(len(spectrum[0]))), np.divide(spectrum[0], np.max(spectrum[0])), lw = 0.5)
plt.xlabel('Energy (keV)')
plt.ylabel('Normalized Counts')
plt.annotate('', xy = (60, 0.5),xytext = (2, 0) ,xycoords = 'data' ,textcoords = 'offset points', arrowprops = dict(arrowstyle = '->'))
plt.annotate('', xy = (59, 0.5), xytext = (2, 0) ,xycoords = 'data' ,textcoords = 'offset points', arrowprops = dict(arrowstyle = '<-'))
plt.text(60.5, 0.52,'FWHM = ' + '%.2f'%(files['Am241'][2] * popt[0] * FWHM))
plt.savefig('/Users/sean/CdTe/gammaFlood/images/20170310_H100_gamma_eu155_-10C_spectrum_keV.0V.eps')
plt.close()


plt.plot(np.multiply(conversion,range(plotBounds[0], plotBounds[1])), np.divide(spectrum[0], np.max(spectrum[0]))[plotBounds[0]:plotBounds[1]], lw = 0.5)
plt.plot(np.multiply(conversion,range(fitBounds[0], fitBounds[1])), np.vectorize(gauss)(range(fitBounds[0], fitBounds[1]), *popt))
plt.xlabel('Energy (keV)')
plt.ylabel('Normalized Counts')
plt.savefig('/Users/sean/CdTe/gammaFlood/images/20170310_H100_gamma_eu155_-10C_fit_HIGH_keV.0V.eps')
plt.close()

### for every event with energy between 40 and 70 keV, add a count at position (rawx, rawy)
pixelMap = np.zeros((32,32))
for i in range(len(events['channel'])):
	if 40 < conversion*events['channel'][i] < 70:
		pixelMap[events['rawx'][i], events['rawy'][i]] += 1

pixelMap[10][10] = 0
pixelMap[0][4] = 0
pixelMap[0][6] = 0

for x in pixelMap:
	print(x)


heatmap = plt.imshow(pixelMap)
plt.colorbar(heatmap, label = 'Counts')
plt.savefig('/Users/sean/CdTe/gammaFlood/images/20170310_H100_gamma_eu155_-10C_heatmap_40-70keV.0V.eps')
plt.close()

plt.hist(pixelMap.flatten(), bins = 30, histtype = 'step')
plt.xlabel('Counts')
plt.ylabel('Pixels')
plt.savefig('/Users/sean/CdTe/gammaFlood/images/20170310_H100_gamma_eu155_-10C_pixelHistogram_40-70keV.0V.eps')
plt.close()'''




