from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import pickle
import skimage.io as imgio
from astropy.modeling import models, fitting

# From http://www.nndc.bnl.gov/nudat2/indx_dec.jsp
Co = 122.06
Am = 59.54
line_high = 10000
line_low = Am*line_high/Co

# Format of the input data
# events = np.array([[{'TIME': [], 'CHANNEL': [], 'GRADE': [], 'STIM': [], 'PH_COM': []} for i in range(32)] for j in range(32)])

# Columns: x, y, gain x13, offset x13 (one for each grade)
rows =[]

gainOffset =  fits.open('/disk/lif2/spike/detectorData/H100/H100_long_gamma_Co57_Am241_-10_gain_offset.0V.fits')[1].data

gain = np.zeros((34, 34))
offset = np.zeros((34, 34))

gain[1:33, 1:33] = gainOffset['GAIN'].reshape(32, 32).T
offset[1:33, 1:33] = gainOffset['OFFSET'].reshape(32, 32).T

'''
print(gainOffset['RAWX'])
print(gainOffset['RAWY'])
print(gainOffset['GAIN'])
print(gain)
'''

for x in range(32):
	for y in range(32):

		row = [x, y, *[1 for i in range(13)], *[0 for j in range(13)]]
		
		dataCo = fits.open('/disk/lif2/spike/detectorData/H100/longGammaFlood/pixelData/H100_long_gamma_Co57_-10_0V_x' + str(x) + '_y' + str(y) + '.fits', memmap=True)[1].data

		dataAm = fits.open('/disk/lif2/spike/detectorData/H100/longGammaFlood/pixelData/H100_long_gamma_Am241_-10_0V_x' + str(x) + '_y' + str(y) + '.fits', memmap=True)[1].data


		if len(dataCo['CHANNEL']) and len(dataAm['CHANNEL']):

			channelGradeCo = [[] for i in range(np.max(dataCo['GRADE'])+1)]
			channelGradeAm = [[] for i in range(np.max(dataAm['GRADE'])+1)]

			for i in range(len(dataCo['CHANNEL'])):
				if (dataCo['CHANNEL'][i] < 20000):
					temp = dataCo['PH_COM'][i].reshape(3,3)
					mask = (temp > 0).astype(int)
					channel = np.sum(np.add(np.multiply(np.multiply(mask, temp), gain[x:x + 3, y:y + 3]), np.multiply(mask, offset[x:x + 3, y:y + 3])))
					channelGradeCo[dataCo['GRADE'][i]].append(channel)

			for i in range(len(dataAm['CHANNEL'])):
				if (dataAm['CHANNEL'][i] < 20000):
					temp = dataAm['PH_COM'][i].reshape(3,3)
					mask = (temp > 0).astype(int)
					channel = np.sum(np.add(np.multiply(np.multiply(mask, temp), gain[x:x + 3, y:y + 3]), np.multiply(mask, offset[x:x + 3, y:y + 3])))
					channelGradeCo[dataAm['GRADE'][i]].append(channel)

			for grade in range(len(channelGradeCo)):

				plt.figure()
				
				if len(channelGradeCo[grade]) and len(channelGradeAm[grade]):

					spectrumCo = np.histogram(channelGradeCo, bins=int(np.ceil(np.max(channelGradeCo))))
					spectrumAm = np.histogram(channelGradeAm, bins=int(np.ceil(np.max(channelGradeAm))))

					plt.plot(range(len(spectrumCo[0])), spectrumCo[0], label = 'Co57')
					plt.plot(range(len(spectrumAm[0])), spectrumAm[0], label = 'Am241', color = 'r')
					plt.ylabel('Counts')
					plt.xlabel('Channel')
					plt.legend()					

					if np.max(channelGradeCo[grade]) > 12000 and np.max(channelGradeAm[grade]) > 6000:
						centroid_high = np.argmax(spectrumCo[0][6000:12000]) + 6000
						fit_channels_high = np.arange(centroid_high-100, centroid_high + 250)
						g_init_high = models.Gaussian1D(amplitude=spectrumCo[0][centroid_high], mean=centroid_high, stddev = 75)
						fit_g_high = fitting.LevMarLSQFitter()
						g_high = fit_g_high(g_init_high, fit_channels_high, spectrumCo[0][fit_channels_high[0]:fit_channels_high[-1]+1])

						centroid_low = np.argmax(spectrumAm[0][4000:6000]) + 4000
						fit_channels_low = np.arange(centroid_low-100, centroid_low + 250)
						g_init_low = models.Gaussian1D(amplitude=spectrumAm[0][centroid_low], mean=centroid_low, stddev = 75)
						fit_g_low = fitting.LevMarLSQFitter()
						g_low = fit_g_low(g_init_low, fit_channels_low, spectrumAm[0][fit_channels_low[0]:fit_channels_low[-1]+1])
						
						row[2 + grade] = (line_high-line_low)/(g_high.mean - g_low.mean)
						row[15 + grade] = (line_low*g_high.mean - line_high*g_low.mean)/(g_high.mean - g_low.mean)
						plt.plot(fit_channels_high, g_high(fit_channels_high))
						plt.plot(fit_channels_low, g_low(fit_channels_low))
						plt.show()

				#plt.savefig('/disk/lif2/spike/detectorData/H100/figures/pixelFits/H100_long_gamma_Co57_Am241_-10_x' + str(x) + '_y' + str(y) + '_gain_offset_grade' + str(grade) + '_linefit.0V.eps')
				plt.close()
		
		rows.append(row)

rows = np.array(rows)
#print(rows)

columns = rows.T
fits_columns = []
col0 = fits.Column(name='RAWX', format='K', array=columns[0])
col1 = fits.Column(name='RAWY', format='K', array=columns[1])
fits_columns.append(col0)
fits_columns.append(col1)
for i in range(13):
	fits_columns.append(fits.Column(name='GAIN_GRADE' + str(i), format='D', array=columns[i + 2]))
for i in range(13):
	fits_columns.append(fits.Column(name='OFFSET_GRADE' + str(i), format='D', array=columns[i + 15]))
t = fits.BinTableHDU.from_columns(fits_columns)
#t.writeto('/Volumes/LaCie/CdTe/longGammaFlood/20170908_H100_long_gamma_Co57_Am241_-10_gain_offset_grade.fits')
