from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
#from matplotlib.axes import Axes


def gauss(x, a, b, c):
	return a*np.exp(np.square(x-b)/((-2.0)*np.square(c)))

def linef(x, m, b):
	return (m*x) + b

def totalf(x, c, d, e, f, g, h, l, m, n, k):
	return (gauss(x, c, d, e) + gauss(x, f, g, h)  + gauss(x, l, m, n) + k)

grv = ['3.35', '3.41', '3.46', '3.51', '3.56', '3.60', '3.66', '3.70', '3.75', '3.80', \
		'3.86', '3.91', '3.96', '4.00', '4.05', '4.10', '4.16']

lineWidth = [[], [], []]

linePos = [[], [], []]

lineStrength = [[], [], []]

widthError = [[], [], []]

posError = [[], [], []]

strengthError = [[], [], []]

for voltage in grv:
	file = fits.open('/Users/sean/Desktop/CdTeData/H100_-400V_-10C_Co57_' + voltage + '.fits')
	#print(file[1].header)
	START = 0
	while file[1].data.field('TEMP')[START] < -50:
		START += 1

	#print(START)
	
	temp = file[1].data.field('PH')[START:]
	file.close()
	nanlist = []
	for i in range(len(temp)):
		if np.isnan(temp[i]) or temp[i]<=0 or temp[i]>1.2E4:
			nanlist.append(i)

	#print(nanlist)
	newtemp = np.delete(temp, nanlist)
	#print(len(newtemp))
	#print(np.max(newtemp))
	#print(newtemp[:20])
	#print(newtemp[0])

	#print(np.min(file[1].data.field('PH')[START:]))
	#print(np.max(file[1].data.field('PH')[START:]))

	spectrum = np.histogram(newtemp, bins = int(np.max(newtemp)))
	
	'''plt.plot(range(len(spectrum[0])), np.divide(spectrum[0], np.max(spectrum[0])), lw = 0.5)
	plt.xlabel('Channel')
	plt.ylabel('Normalized Counts')
	plt.savefig('/Users/sean/Desktop/CdTeData/images/H100_-400V_-10C_Co57_' + voltage + '.eps')
	plt.close()'''
	
	popt, pcov = curve_fit(totalf, range(8500, 10500), np.divide(spectrum[0], np.max(spectrum[0]))[8500:10500], p0 = [1.0, 9250, 1.0, 1.0, 9750, 1.0, 1.0, 9600, 1.0, 0.0])

	# stdev error
	errors = np.sqrt(np.diag(pcov))

	for i in range(len(lineStrength)):
		lineStrength[i].append(popt[3*i])
		linePos[i].append(popt[(3*i) + 1])
		lineWidth[i].append(popt[(3*i) + 2])
		strengthError[i].append(errors[3*i])
		posError[i].append(errors[(3*i) + 1])
		widthError[i].append(errors[(3*i) + 2])

	'''plt.plot(range(8500, 10500), np.divide(spectrum[0], np.max(spectrum[0]))[8500:10500], lw = 0.5)
	plt.plot(range(8500, 10500), np.vectorize(totalf)(range(8500, 10500), *popt))
	plt.xlabel('Channel')
	plt.ylabel('Normalized Counts')
	plt.savefig('/Users/sean/Desktop/CdTeData/images/H100_-400V_-10C_Co57_' + voltage + 'LineFit.eps')
	plt.close()
'''
colors = ['blue', 'orange', 'green']
for i in range(len(linePos)):
	plt.figure()
	plt.errorbar(np.vectorize(float)(grv), linePos[i], yerr = posError[i], color = colors[i])
	plt.xlabel('Guard Ring Voltage (V)')
	plt.ylabel('Line Position (Channel)')
	plt.savefig('/Users/sean/Desktop/CdTeData/images/H100_-400V_-10C_Co57_LinePos_' + str(int(np.mean(linePos[i]))) + '.eps')
	plt.close()
'''
plt.figure()
plt.errorbar(np.vectorize(float)(grv), np.sum(lineWidth, axis = 0), yerr = np.sqrt(np.sum(np.square(widthError), axis = 0)))
plt.xlabel('Guard Ring Voltage (V)')
plt.ylabel('Line Width Sum (Channels)')
plt.savefig('/Users/sean/Desktop/CdTeData/images/H100_-400V_-10C_Co57_LineWidths_Sum.eps')
plt.close()


plt.figure()
for i in range(len(linePos)):
	plt.errorbar(np.vectorize(float)(grv), linePos[i], yerr = posError[i])'''

