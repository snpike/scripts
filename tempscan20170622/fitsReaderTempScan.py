from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.io import readsav
import pickle

highELine = 59.5412

# Gain and offset from fitting three different lines at -10C
channelToEnergy =[ 0.03987333,  0.34386187]

def gauss(x, a, b, c):
	pdf =  a*np.exp(np.square(x-b)/((-2.0)*np.square(c)))
	return pdf

lineWidth = []

linePos = []

widthError = []

posError = []
fitBounds = [1475, 1550]
plotBounds = [1300, 1600]
FWHM = 2.0*np.sqrt(np.log(4.0))

'''
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

#print(gain)
#print(offset)

file = fits.open('/Users/sean/CdTe/tempscan20170622/20170622_H100_slowtempscan_Am241.0V.fits')
data = file[1].data
file.close()
START = 0
i = 0

# Skip the beginning
while data.field('TEMP')[i] < -50:
	i += 1
START = i

events = [{'channel': [], 'rawx': [], 'rawy': [], 'grade': [], 'temperature': [], 'time': []} for j in range(9)]

j = 0
print(j)
while i < len(data.field('PH')):
	rawx = data.field('RAWX')[i]
	rawy = data.field('RAWY')[i]
	grade = data.field('GRADE')[i]
	T = data.field('TEMP')[i]
	time = data.field('TIME')[i]
	if T < (-2*j)-1:
		j+=1
		print(j)
	if ((-2*j) + 0.4) < T < ((-2*j)+0.6):
		if grade == 0:
		### Apply the gain and offset ONLY to pixels with positive counts
			temp = data.field('PH_COM')[i].reshape(3,3)
			if np.sum(temp) > 0:
				mask = (temp > 0).astype(int)
				channel = np.sum(np.add(np.multiply(np.multiply(mask, temp), gain[rawx:rawx + 3, rawy:rawy + 3]), np.multiply(mask, offset[rawx:rawx + 3, rawy:rawy + 3])))
				#if grade < 5:
				channel = (channel * PHGainOffset[rawx][rawy]['gain'][grade]) + PHGainOffset[rawx][rawy]['offset'][grade]
				if (not np.isnan(channel)) and (channel < 3500):
					events[j]['channel'].append(channel)
					events[j]['rawx'].append(rawx)
					events[j]['rawy'].append(rawy)
					events[j]['grade'].append(grade)
					events[j]['time'].append(time)
					events[j]['temperature'].append(T)
	i += 1

for x in events:
	dumpFile = open('/Users/sean/CdTe/tempscan20170622/20170622_H100_slowtempscan_Am241_' + '%.1f'%(np.mean(x['temperature'])) + 'SingleEvent_dump.0V.pkl', 'wb')
	pickle.dump(x, dumpFile)
	dumpFile.close()

'''
tempLabels = ['0.5', '-1.5', '-3.5', '-5.5', '-7.5', '-9.5', '-11.5', '-13.5', '-15.5']

events = []
for x in tempLabels:
	pFile = open('/Users/sean/CdTe/tempscan20170622/20170622_H100_slowtempscan_Am241_' + x + 'SingleEvent_dump.0V.pkl', 'rb')
	events.append(pickle.load(pFile))
	pFile.close()

conversionList = []
#error in channel -> energy gain
conversionErrList = []
meanTemp = []
tempErr = []
for x in events:
	meanTemp.append(np.mean(x['temperature']))
	tempErr.append(np.std(x['temperature']))
	spectrum = np.histogram(x['channel'], bins = int(np.max(x['channel'])))
	
	plt.plot(range(len(spectrum[0])), np.divide(spectrum[0], np.max(spectrum[0])), lw = 1.0)
	plt.xlabel('Channel')
	plt.ylabel('Normalized Counts')
	plt.savefig('/Users/sean/CdTe/tempscan20170622/images/20170622_H100_slowtempscan_Am241_' + '%.1f'%(meanTemp[-1]) + '_SingleEvent_spectrum.0V.eps')
	plt.close()
	centroid = np.argmax(spectrum[0])
	fitBounds = [centroid - 7, centroid + 15]
	popt, pcov = curve_fit(gauss, range(fitBounds[0], fitBounds[1]), np.divide(spectrum[0], np.max(spectrum[0]))[fitBounds[0]:fitBounds[1]], p0 = [1.0, centroid, 1.0])
	conversion = highELine/popt[1]
	
	#stdev error
	errors = np.sqrt(np.diag(pcov))
	lineWidth.append(conversion*FWHM*popt[2])

	conversionErr = errors[1]*highELine/np.square(popt[1])

	conversionList.append(conversion)
	conversionErrList.append(conversionErr)
	widthError.append(FWHM*np.sqrt(np.square(errors[2]*conversion) + np.square(conversionErr*popt[2])))
	
	plt.plot(np.multiply(conversion,range(len(spectrum[0]))), np.divide(spectrum[0], np.max(spectrum[0])), lw = 1.0)
	plt.xlabel('Energy (keV)')
	plt.ylabel('Normalized Counts')
	plt.savefig('/Users/sean/CdTe/tempscan20170622/images/20170622_H100_slowtempscan_Am241_' + '%.1f'%(meanTemp[-1]) + '_SingleEvent_spectrum_keV.0V.eps')
	plt.close()

	
	plt.plot(np.multiply(conversion,range(plotBounds[0], plotBounds[1])), np.divide(spectrum[0], np.max(spectrum[0]))[plotBounds[0]:plotBounds[1]], lw = 1.0)
	plt.plot(np.multiply(conversion,range(fitBounds[0], fitBounds[1])), np.vectorize(gauss)(range(fitBounds[0], fitBounds[1]), *popt))
	plt.xlabel('Energy (keV)')
	plt.ylabel('Normalized Counts')
	plt.savefig('/Users/sean/CdTe/tempscan20170622/images/20170622_H100_slowtempscan_Am241_' + '%.1f'%(meanTemp[-1]) + '_SingleEvent_fit_HIGH_keV.0V.eps')
	plt.close()
	
	'''
	# fit the trigger signal to a gaussian
	centroid = np.argmax(spectrum[0][2300:2380]) + 2300
	fitBounds = [centroid-20, centroid+20]
	plotBounds = [centroid-40, centroid+40]
	popt, pcov = curve_fit(gauss, range(fitBounds[0], fitBounds[1]), np.divide(spectrum[0], np.max(spectrum[0]))[fitBounds[0]:fitBounds[1]], p0 = [1.0, centroid, 1.0])

	lineWidth.append(conversion*FWHM*popt[2])
	widthError.append(FWHM*np.sqrt(np.square(errors[2]*conversion) + np.square(conversionErr*popt[2])))

	plt.plot(np.multiply(conversion,range(plotBounds[0], plotBounds[1])), np.divide(spectrum[0], np.max(spectrum[0]))[plotBounds[0]:plotBounds[1]], lw = 0.5)
	plt.plot(np.multiply(conversion,range(fitBounds[0], fitBounds[1])), np.vectorize(gauss)(range(fitBounds[0], fitBounds[1]), *popt))
	plt.xlabel('Energy (keV)')
	plt.ylabel('Normalized Counts')
	plt.savefig('/Users/sean/CdTe/tempscan20170622/images/20170622_H100_slowtempscan_Am241_' + '%.1f'%(meanTemp[-1]) + '_fit_TRIGGER_keV.0V.eps')
	plt.close()'''
	

plt.figure()
plt.errorbar(meanTemp, lineWidth, yerr = widthError, xerr = tempErr, fmt = '.', elinewidth = 1.0)
plt.xlabel('Temperature (C)')
plt.ylabel('FWHM (keV)')
plt.savefig('/Users/sean/CdTe/tempscan20170622/images/lineWidths_SingleEvent_keV.0V.eps')
plt.close()


plt.figure()
plt.errorbar(meanTemp, conversionList, yerr = conversionErrList, xerr = tempErr, fmt = '.', elinewidth = 1.0)
#plt.errorbar(-9.721488, channelToEnergy[0], fmt = '.')
plt.xlabel('Temperature (C)')
plt.ylabel('Channel to keV gain (kev/channel)')
plt.savefig('/Users/sean/CdTe/tempscan20170622/images/channel2keV_SingleEvent.0V.eps')
plt.close()
'''
meanTemp = []
tempErr = []
gradeHist = [np.histogram(events[j]['grade'], bins = 13)[0] for j in range(9)]
width = 0.5

pixelNum = [[] for i in range(9)]

for j in range(9):
	pixelNum[j].append(gradeHist[j][0])
	pixelNum[j].append(np.sum(gradeHist[j][1:4]))
	pixelNum[j].append(np.sum(gradeHist[j][5:8]))
	pixelNum[j].append(np.sum(gradeHist[j][9:12]))
	meanTemp.append(np.mean(events[j]['temperature']))
	tempErr.append(np.std(events[j]['temperature']))

numError = [[] for i in range(9)]
for j in range(9):
	for i in range(len(pixelNum[j])):
		numError[j].append(100.0*np.sqrt(1.0 - (pixelNum[j][i]/np.sum(pixelNum[j])))*np.sqrt(pixelNum[j][i])/np.sum(pixelNum[j]))

for j in range(9):
	pixelNum[j] = np.divide(pixelNum[j], np.sum(np.divide(pixelNum[j], 100)))


#print(pixelNum)

for j in range(9):
	fig, ax = plt.subplots()
	ax.set_xticks(range(4))
	ax.set_ylabel('Percentage of Events')
	ax.set_xlabel('Number of Pixels')
	rects = ax.bar(range(4), pixelNum[j], width, color='r')
	ax.set_xticklabels(('1', '2', '3', '4'))
	plt.savefig('/Users/sean/CdTe/tempscan20170622/images/20170622_H100_slowtempscan_Am241_' + '%.1f'%(meanTemp[j]) + '_PixelNum.eps')
	plt.close()

colors = ['r','b','g','m']
f, axes = plt.subplots(4, sharex=True, sharey=False)
for i in range(4):
	axes[i].errorbar(meanTemp, np.transpose(pixelNum)[i], xerr = tempErr, yerr = np.transpose(numError)[i], lw = 1.0, color = colors[i])
	axes[i].set_ylim([np.mean(np.transpose(pixelNum)[i]) - 0.75, np.mean(np.transpose(pixelNum)[i]) + 0.75])
axes[-1].set_xlabel('Temperature (C)')
axes[0].set_ylabel('% 1-Pixel')
axes[1].set_ylabel('% 2-Pixel')
axes[2].set_ylabel('% 3-Pixel')
axes[3].set_ylabel('% 4-Pixel')
f.subplots_adjust(hspace=0)
plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
plt.savefig('/Users/sean/CdTe/tempscan20170622/images/20170622_H100_slowtempscan_Am241_EventPercents.eps')
plt.close()

heatmap = plt.imshow(np.histogram2d(events[5]['rawx'], events[5]['rawy'], bins = [32,32])[0])
plt.colorbar(heatmap, label = 'Counts')
plt.savefig('/Users/sean/CdTe/tempscan20170622/images/20170622_H100_slowtempscan_Am241_' + '%.1f'%(meanTemp[5]) + 'pixelMap.0V.eps')
plt.close()'''


