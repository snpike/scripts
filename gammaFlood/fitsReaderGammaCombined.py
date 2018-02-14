from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.io import readsav
from scipy.special import erf
import pickle


def gauss(x, a, b, c):
	pdf =  a*np.exp(np.square(x-b)/((-2.0)*np.square(c)))
	#cdf = (1.0 + erf(alpha*(x-b)/(np.sqrt(2.0)*c)))/2.0
	#return pdf*cdf + k
	return pdf

def line(x, m, b):
	return (m*x) + b

FWHM = 2.0*np.sqrt(np.log(4.0))

files = {'Am241': ['20170315_H100_gamma_Am241_', 59.5412], 'Co57': ['20170313_H100_gamma_Co57_', 122.0614], 'Eu155': ['20170310_H100_gamma_eu155_', 86.545]}
centroids = []
centroidErr = []
for x in files:
	events = pickle.load(open('/Users/sean/CdTe/tempscan20170622/gammaFlood/' + files[x][0] + '-10C_dump.0V.pkl', 'rb'))

	spectrum = np.histogram(events['channel'], bins = int(np.max(events['channel'])))

	#plt.plot(range(len(spectrum[0])), np.divide(spectrum[0], np.max(spectrum[0])), lw = 0.5)

	centroid = np.argmax(spectrum[0])
	fitBounds = [centroid - 7, centroid + 15]
	plotBounds = [fitBounds[0] - 50, fitBounds[1] + 42]
	popt, pcov = curve_fit(gauss, range(fitBounds[0], fitBounds[1]), np.divide(spectrum[0], np.max(spectrum[0]))[fitBounds[0]:fitBounds[1]], p0 = [1.0, centroid, 1.0])
	centroids.append(popt[1])
	errors = np.sqrt(np.diag(pcov))
	centroidErr.append(errors[1])
	files[x].append(popt[2])
	files[x].append(errors[2])

	#stdev error
	#errors = np.sqrt(np.diag(pcov))

	#conversionErr = errors[1]*highELine/np.square(popt[1])

	#widthError = FWHM*np.sqrt(np.square(errors[2]*conversion) + np.square(conversionErr*popt[2]))


popt, pcov = curve_fit(line, centroids, [files[x][1] for x in files])

for x in files:
	events = pickle.load(open('/Users/sean/CdTe/tempscan20170622/gammaFlood/' + files[x][0] + '-10C_dump.0V.pkl', 'rb'))

	spectrum = np.histogram(events['channel'], bins = int(np.max(events['channel'])))

	plt.plot(np.vectorize(line)(range(len(spectrum[0])), *popt), np.divide(spectrum[0], np.max(spectrum[0])), lw = 1.0, label = x)
	#plt.legend(label = x)

plt.xlabel('Energy (keV)')
plt.ylabel('Normalized Counts')
plt.annotate('', xy = (60, 0.5),xytext = (2, 0) ,xycoords = 'data' ,textcoords = 'offset points', arrowprops = dict(arrowstyle = '->'))
plt.annotate('', xy = (59, 0.5), xytext = (2, 0) ,xycoords = 'data' ,textcoords = 'offset points', arrowprops = dict(arrowstyle = '<-'))
plt.text(28, 0.52,'FWHM = ' + '%.2f'%(files['Am241'][2] * popt[0] * FWHM))
plt.annotate('', xy = (87, 0.5),xytext = (2, 0) ,xycoords = 'data' ,textcoords = 'offset points', arrowprops = dict(arrowstyle = '->'))
plt.annotate('', xy = (86, 0.5), xytext = (2, 0) ,xycoords = 'data' ,textcoords = 'offset points', arrowprops = dict(arrowstyle = '<-'))
plt.text(87.5, 0.52,'%.2f'%(files['Eu155'][2] * popt[0] * FWHM))
plt.annotate('', xy = (122.5, 0.5),xytext = (2, 0) ,xycoords = 'data' ,textcoords = 'offset points', arrowprops = dict(arrowstyle = '->'))
plt.annotate('', xy = (121.5, 0.5), xytext = (2, 0) ,xycoords = 'data' ,textcoords = 'offset points', arrowprops = dict(arrowstyle = '<-'))
plt.text(122.7, 0.52,'%.1f'%(files['Co57'][2] * popt[0] * FWHM))
plt.legend()
plt.savefig('/Users/sean/CdTe/tempscan20170622/images/Am241_Co57_Eu155_-10C_combinedSpectrum_keV.0V.eps')
plt.close()

print(files['Am241'][2] * popt[0] * FWHM)
print(files['Eu155'][2] * popt[0] * FWHM)
print(files['Co57'][2] * popt[0] * FWHM)

events = pickle.load(open('/Users/sean/CdTe/tempscan20170622/gammaFlood/' + files['Am241'][0] + '-10C_dump.0V.pkl', 'rb'))
pixelMap = np.zeros((32,32))
for i in range(len(events['channel'])):
	if 40 < line(events['channel'][i], *popt) < 70:
		pixelMap[events['rawx'][i], events['rawy'][i]] += 1

pixelMap[10][10] = 0
pixelMap[0][4] = 0
pixelMap[0][6] = 0

heatmap = plt.imshow(pixelMap)
plt.colorbar(heatmap, label = 'Counts')
plt.savefig('/Users/sean/CdTe/tempscan20170622/images/20170315_H100_gamma_Am241_-10C_heatmap.0V.eps')
plt.close()

plt.hist(pixelMap.flatten(), bins = 30, histtype = 'step')
plt.xlabel('Counts')
plt.ylabel('Pixels')
plt.savefig('/Users/sean/CdTe/tempscan20170622/images/20170315_H100_gamma_Am241_-10C_pixelHistogram_40-70keV.0V.eps')
plt.close()
