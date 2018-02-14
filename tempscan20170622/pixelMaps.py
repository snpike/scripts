from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.io import readsav
from scipy.special import erf
import pickle
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

file = fits.open('/Users/sean/CdTe/tempscan20170622/20170622_H100_slowtempscan_Am241.fits')
i = 0
#PH = []
# Skip the beginning
while file[1].data.field('TEMP')[i] < -50:
	i += 1
temps = [[] for j in range(9)]
#pixelMaps = [np.zeros((32,32)) for j in range(9)]
exposureT = [0.0 for j in range(9)]

j = 0
print(j)
while i < len(file[1].data.field('PH')) - 1:
	#rawx = file[1].data.field('RAWX')[i]
	#rawy = file[1].data.field('RAWY')[i]
	#grade = file[1].data.field('GRADE')[i]
	T = file[1].data.field('TEMP')[i]
	if T < (-2*j)-1:
		j+=1
		print(j)
	if ((-2*j) + 0.4) < T < ((-2*j)+0.6):
		if not (((-2*j) + 0.4) < file[1].data.field('TEMP')[i-1] < ((-2*j)+0.6)):
			START = file[1].data.field('TIME')[i]

		if not (((-2*j) + 0.4) < file[1].data.field('TEMP')[i+1] < ((-2*j)+0.6)):
			exposureT[j] += file[1].data.field('TIME')[i] - START

		### Apply the gain and offset ONLY to pixels with positive counts
		temp = file[1].data.field('PH_COM')[i].reshape(3,3)
		if np.sum(temp) > 0:
			mask = (temp > 0).astype(int)
			energy = np.sum(np.multiply(np.add(np.multiply(mask, offset[rawx:rawx + 3, rawy:rawy + 3]), temp), gain[rawx:rawx + 3, rawy:rawy + 3]))
			if grade < 5:
				energy = (energy + PHGainOffset[rawx][rawy]['offset'][grade]) * PHGainOffset[rawx][rawy]['gain'][grade]
			if (not np.isnan(energy)) and (energy > 10):
				pixelMaps[j][rawx,rawy] += 1
				temps[j].append(T)
	i += 1
exposureT[j] += file[1].data.field('TIME')[i] - START

file.close()
print(exposureT)
#pickle.dump(pixelMaps, open('/Users/sean/CdTe/tempscan20170622/20170622_H100_slowtempscan_Am241_PixelMaps.pkl', 'wb'))
#pickle.dump(temps, open('/Users/sean/CdTe/tempscan20170622/20170622_H100_slowtempscan_Am241_PixelMapTemps.pkl', 'wb'))
'''
pixelMaps = pickle.load(open('/Users/sean/CdTe/tempscan20170622/20170622_H100_slowtempscan_Am241_PixelMaps.pkl', 'rb'))
temps = pickle.load(open('/Users/sean/CdTe/tempscan20170622/20170622_H100_slowtempscan_Am241_PixelMapTemps.pkl', 'rb'))

for i in range(len(pixelMaps)):
	pixelMaps[i][10,10] = 0

exposureT = [7368.0, 6422.0, 6442.0, 6436.0, 6559.0, 6702.0, 6708.0, 6627.0, 6555.0]

# Plotting and fitting from here
meanTemp = []
tempErr = []
for j in range(9):
	meanTemp.append(np.mean(temps[j]))
	tempErr.append(np.std(temps[j]))
	
	heatmap = plt.imshow(np.divide(pixelMaps[j], exposureT[j]))
	plt.colorbar(heatmap, label = 'Counts per second')
	plt.savefig('/Users/sean/CdTe/tempscan20170622/images/20170622_H100_slowtempscan_Am241_' + '%.1f'%(meanTemp[j]) + '_heatmap.eps')
	plt.close()

plt.plot(meanTemp, np.divide([np.sum(x) for x in pixelMaps], exposureT))
plt.show()
plt.close()

