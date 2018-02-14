from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.io import readsav
import pickle
from scipy import interpolate
from matplotlib.animation import MovieWriter

highELine = 122.0614

def gauss(x, a, b, c):
	pdf =  a*np.exp(np.square(x-b)/((-2.0)*np.square(c)))
	#return pdf*cdf + k
	return pdf

FWHM = 2.0*np.sqrt(np.log(4.0))

countRate = 0.18

### Load the pickled data and do some plotting/fit the 59 keV line to a gaussian. Use this fit to convert from channels to energy.

events = pickle.load(open('/Users/sean/CdTe/gammaFlood/20170313_H100_gamma_Co57_-10C_dump.0V.pkl', 'rb'))

spectrum = np.histogram(events['channel'], bins = range(3501))

simPDF = np.divide(spectrum[0], np.sum(spectrum[0]))
simCDF = [np.sum(simPDF[:i]) for i in range(len(simPDF))]
simCDF[-1] = 1.0
inverseCDF = interpolate.interp1d(simCDF, range(len(simCDF)))

n = 100
simEvents = []
fitError = []
counts = []
for j in range(n):
	simEvents.append(inverseCDF(np.random.random(1000)))
	simSpec = np.histogram(simEvents, bins = range(3501))[0]
	centroid = 3054
	fitBounds = [centroid - 7, centroid + 15]
	popt, pcov = curve_fit(gauss, range(fitBounds[0], fitBounds[1]), simSpec[fitBounds[0]:fitBounds[1]], p0 = [1.0, centroid, 1.0])

	#stdev error
	errors = np.sqrt(np.diag(pcov))
	counts.append(len(simEvents) * 1000)
	fitError.append(errors[1])
	plt.figure()
	plt.xlabel('Channel')
	plt.ylabel('Counts')
	plt.plot(range(3500), simSpec)
	plt.savefig('/Users/sean/CdTe/gammaFlood/images/movieFrames/gammaSimCo_step' + str(j+1) + '.png')
	plt.close()
	if ((j+1) % 10 == 0):
		print(popt[1])
		plt.figure()
		plt.xlabel('Channel')
		plt.ylabel('Counts')
		plt.plot(range(3500), simSpec)
		plt.savefig('/Users/sean/CdTe/gammaFlood/images/gammaSimCo_step' + str(j+1) + '.eps')
		plt.close()
plt.figure()
plt.xlabel('Total Counts')
plt.ylabel('Centroid Error')
plt.plot(counts[5:], fitError[5:])
plt.savefig('/Users/sean/CdTe/gammaFlood/images/gammaSimCo_error.eps')
plt.close()

