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


gradeHist = [[0 for j in range(13)] for k in range(9)]
file = fits.open('/Users/sean/CdTe/tempscan20170622/20170622_H100_slowtempscan_Am241.fits')
data = file[1].data
file.close()
i = 0

# Skip the beginning
while data.field('TEMP')[i] < -50:
	i += 1

j = 0
print(j)
while i < len(data.field('PH')) - 1:
	rawx = data.field('RAWX')[i]
	rawy = data.field('RAWY')[i]
	grade = data.field('GRADE')[i]
	T = data.field('TEMP')[i]
	if T < (-2*j)-1:
		j+=1
		print(j)
	if ((-2*j) + 0.4) < T < ((-2*j)+0.6):
		### Apply the gain and offset ONLY to pixels with positive counts
		temp = data.field('PH_COM')[i].reshape(3,3)
		if np.sum(temp) > 0:
			mask = (temp > 0).astype(int)
			energy = np.sum(np.multiply(np.add(np.multiply(mask, offset[rawx:rawx + 3, rawy:rawy + 3]), temp), gain[rawx:rawx + 3, rawy:rawy + 3]))
			if grade < 5:
				energy = (energy + PHGainOffset[rawx][rawy]['offset'][grade]) * PHGainOffset[rawx][rawy]['gain'][grade]
			if (not np.isnan(energy)) and (energy > 10):
				gradeHist[j][grade] += 1
				
	i += 1
'''
temps = pickle.load(open('/Users/sean/CdTe/tempscan20170622/20170622_H100_slowtempscan_Am241_PixelMapTemps.pkl', 'rb'))

exposureT = [7368.0, 6422.0, 6442.0, 6436.0, 6559.0, 6702.0, 6708.0, 6627.0, 6555.0]

# Plotting and fitting from here
meanTemp = []
tempErr = []
gradeHist = [[1310416, 132397, 131471, 129318, 131182, 6266, 5684, 5468, 5403, 16888, 15236, 15416, 15934], [1129041, 113615, 113954, 111156, 113927, 5770, 5143, 5067, 5029, 16139, 14863, 14719, 15417], [1128182, 112519, 114230, 111345, 113334, 5937, 5305, 5239, 5017, 17243, 15680, 15565, 16153], [1125684, 111101, 112594, 109556, 112315, 6013, 5387, 5420, 5049, 17853, 16507, 16441, 17229], [1147255, 112025, 113679, 110216, 113542, 6258, 5565, 5538, 5261, 18915, 17358, 17429, 18257], [1169737, 113983, 114848, 111622, 115432, 6580, 5632, 5803, 5605, 20332, 18570, 18641, 19428], [1168599, 112842, 113878, 110652, 114330, 6837, 5888, 5875, 5642, 20916, 19061, 19086, 20024], [1153432, 110541, 112121, 107669, 111837, 6821, 5908, 5788, 5444, 21241, 19379, 19472, 20241], [1140174, 108208, 109801, 105965, 109821, 6882, 5916, 5885, 5538, 21376, 19622, 19167, 20470]]
width = 0.5

pixelNum = [[] for i in range(9)]

for j in range(9):
	pixelNum[j].append(gradeHist[j][0])
	pixelNum[j].append(np.sum(gradeHist[j][1:4]))
	pixelNum[j].append(np.sum(gradeHist[j][5:8]))
	pixelNum[j].append(np.sum(gradeHist[j][9:12]))
	meanTemp.append(np.mean(temps[j]))
	tempErr.append(np.std(temps[j]))

numError = [[] for i in range(9)]
for j in range(9):
	for i in range(len(pixelNum[j])):
		numError[j].append(100.0*np.sqrt(1.0 - (pixelNum[j][i]/np.sum(pixelNum[j])))*np.sqrt(pixelNum[j][i])/np.sum(pixelNum[j]))

for j in range(9):
	pixelNum[j] = np.divide(pixelNum[j], np.sum(np.divide(pixelNum[j], 100)))


#print(pixelNum)
'''
for j in range(9):
	fig, ax = plt.subplots()
	ax.set_xticks(range(4))
	ax.set_ylabel('Percentage of Events')
	ax.set_xlabel('Number of Pixels')
	rects = ax.bar(range(4), pixelNum[j], width, color='r')
	ax.set_xticklabels(('1', '2', '3', '4'))
	plt.savefig('/Users/sean/CdTe/tempscan20170622/images/20170622_H100_slowtempscan_Am241_' + '%.1f'%(meanTemp[j]) + '_PixelNum.eps')
	plt.close()'''

colors = ['r','b','g','m']
f, axes = plt.subplots(4, sharex=True, sharey=False)
for i in range(4):
	axes[i].errorbar(meanTemp, np.transpose(pixelNum)[i], xerr = tempErr, yerr = np.transpose(numError)[i], lw = 0.5, color = colors[i])
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

