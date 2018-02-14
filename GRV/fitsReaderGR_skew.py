from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.io import readsav
from scipy.special import erf
import pickle
from scipy.stats import invgauss

'''
def gauss(x, a, b, mu, lmbda):
	x = np.abs(a - x) * b
	#boolx = (x > 0).astype(np.float)
	pdf = np.sqrt(lmbda/(2.0 * np.pi * np.power(x, 3))) * np.exp(-lmbda*np.square(x-mu)/(2.0*np.square(mu)*x))
	return pdf'''

def gauss(x, a, b, mu, c, k):
	x = (a - x) * b
	boolx = (x > 0)
	pdf = boolx * np.multiply(np.divide(mu, np.sqrt(2.0 * np.pi * np.power(x, 3))), np.exp(np.square(x-mu)/(-2.0*x)))
	return np.add(np.multiply(c, pdf), k)

grv = ['3.35', '3.41', '3.46', '3.51', '3.56', '3.60', '3.66', '3.70', '3.75', '3.80', \
		'3.86', '3.91', '3.96', '4.00', '4.05', '4.10', '4.16']

lineWidth = []

linePos = []

widthError = []

posError = []

gainOffset =  readsav('/Users/sean/Desktop/CdTeData/h100_-10C_500V_gain.sav')['calstruct']
PHGainOffset = readsav('/Users/sean/Desktop/CdTeData/h100_-10C_500V_fit_2.sav')['calstruct3'].reshape(32,32)

tempgain = []
tempoffset = []

fitBounds = [2500, 3100]
plotBounds = [2500, 3100]
centroid = 3040


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

#print(gain)
#print(offset)



for voltage in grv:
	'''
	file = fits.open('/Users/sean/Desktop/CdTeData/GRV/GRV_fits/H100_-400V_-10C_Co57_' + voltage + '.fits')
	START = 0
	i = 0
	PH = []
	# Skip the beginning
	while file[1].data.field('TEMP')[i] < -50:
		i += 1
	START = i
	# Apply gain/offset correction to PREPHAS and POSTPHAS data. Subtract to get event
	while i < len(file[1].data.field('PH')):
		rawx = file[1].data.field('RAWX')[i]
		rawy = file[1].data.field('RAWY')[i]
		grade = file[1].data.field('GRADE')[i]

		### Apply the gain and offset ONLY to pixels with positive counts
		temp = file[1].data.field('PH_COM')[i].reshape(3,3)
		if np.sum(temp) > 0:
			mask = (temp > 0).astype(int)
			PH.append(np.sum(np.multiply(np.add(np.multiply(mask, offset[rawx:rawx + 3, rawy:rawy + 3]), temp), gain[rawx:rawx + 3, rawy:rawy + 3])))
			if grade < 5:
				PH[-1] = (PH[-1] + PHGainOffset[rawx][rawy]['offset'][grade]) * PHGainOffset[rawx][rawy]['gain'][grade]
		i += 1

	file.close()
	nanlist = []
	for i in range(len(PH)):
		if np.isnan(PH[i]) or PH[i]< 10 or PH[i]>3750:
			nanlist.append(i)

	#print(nanlist)
	newPH = np.delete(PH, nanlist)

	# Plotting and fitting from here
	spectrum = np.histogram(newPH, bins = int(np.max(newPH)))
	pickle.dump(spectrum, open('/Users/sean/Desktop/CdTeData/GRV/GRV_pkl/H100_-400V_-10C_Co57_' + voltage + '_GainCorrectedSpectrum.pkl', 'wb'))
	'''
	pFile = open('/Users/sean/Desktop/CdTeData/GRV/GRV_pkl/H100_-400V_-10C_Co57_' + voltage + '_GainCorrectedSpectrum.pkl', 'rb')
	spectrum = pickle.load(pFile)
	pFile.close()
	'''
	plt.plot(range(len(spectrum[0])), np.divide(spectrum[0], np.max(spectrum[0])), lw = 0.5)
	plt.xlabel('Channel')
	plt.ylabel('Normalized Counts')
	#plt.show()
	plt.savefig('/Users/sean/Desktop/CdTeData/GRV/images/gain_corrected_spectra/H100_-400V_-10C_Co57_' + voltage + '_GAIN.eps')
	plt.close()'''
	gaussOffset = np.mean(np.divide(spectrum[0], np.max(spectrum[0]))[2600:2700])
	popt, pcov = curve_fit(gauss, range(fitBounds[0], fitBounds[1]), np.divide(spectrum[0], np.max(spectrum[0]))[fitBounds[0]:fitBounds[1]], p0 = [3200, 0.005, 0.01, 3.0, gaussOffset], bounds = ([3000, 0.001, 0.01, 0.01, gaussOffset-0.001],[3200, 0.1, 10.0, 10.0, gaussOffset + 0.001]))
	
	# stdev error
	errors = np.sqrt(np.diag(pcov))
	
	print(popt)
	print(pcov)
	print('\n')

	lineWidth.append(popt[0])
	linePos.append(popt[1])
	widthError.append(errors[0])
	posError.append(errors[1])

	
	plt.plot(range(plotBounds[0], plotBounds[1]), np.divide(spectrum[0], np.max(spectrum[0]))[plotBounds[0]:plotBounds[1]], lw = 0.5)
	plt.plot(range(plotBounds[0], plotBounds[1]), np.vectorize(gauss)(range(plotBounds[0], plotBounds[1]), *popt))
	plt.xlabel('Channel')
	plt.ylabel('Normalized Counts')
	#plt.show()
	plt.savefig('/Users/sean/Desktop/CdTeData/GRV/images/LineFit_corrected/H100_-400V_-10C_Co57_' + voltage + 'LineFit_SKEW.eps')
	plt.close()
	


plt.figure()
plt.errorbar(np.vectorize(float)(grv), lineWidth, yerr = widthError)
plt.xlabel('Guard Ring Voltage (V)')
plt.ylabel('Line Width (Channels)')
plt.savefig('/Users/sean/Desktop/CdTeData/GRV/images/H100_-400V_-10C_Co57_LineWidths_SKEW.eps')
plt.close()

plt.figure()
plt.errorbar(np.vectorize(float)(grv), linePos, yerr = posError)
plt.xlabel('Guard Ring Voltage (V)')
plt.ylabel('Line Position (Channel)')
plt.savefig('/Users/sean/Desktop/CdTeData/GRV/images/H100_-400V_-10C_Co57_LinePos_SKEW.eps')
plt.close()
