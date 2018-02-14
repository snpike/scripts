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

for x in range(32):
	for y in range(32):
		row = [x, y, *[1 for i in range(13)], *[0 for j in range(13)]]
		#print(row)
		data = pickle.load(open('/Volumes/LaCie/CdTe/longGammaFlood/pixelData/H100_long_gamma_Co57_Am241_-10_dump_x' + str(x) + '_y' + str(y) + '.0V.pkl', 'rb'))
		if len(data['CHANNEL']):
			channelGrade = [[] for i in range(np.max(data['GRADE'])+1)]
			#print(np.max(data['GRADE']))
			for i in range(len(data['CHANNEL'])):
				if (data['CHANNEL'][i] < 20000):
					channelGrade[data['GRADE'][i]].append(data['CHANNEL'][i])
			for grade in range(len(channelGrade)):
				if len(channelGrade[grade]):
					spectrum = np.histogram(channelGrade[grade], bins=int(np.ceil(np.max(channelGrade[grade]))))
					if np.max(channelGrade[grade]) > 12000:
						centroid_high = np.argmax(spectrum[0][6000:12000]) + 6000
						fit_channels_high = np.arange(centroid_high-100, centroid_high + 250)
						g_init_high = models.Gaussian1D(amplitude=spectrum[0][centroid_high], mean=centroid_high, stddev = 75)
						fit_g_high = fitting.LevMarLSQFitter()
						g_high = fit_g_high(g_init_high, fit_channels_high, spectrum[0][fit_channels_high[0]:fit_channels_high[-1]+1])

						centroid_low = np.argmax(spectrum[0][101:650]) + 101
						fit_channels_low = np.arange(centroid_low-100, centroid_low + 250)
						g_init_low = models.Gaussian1D(amplitude=spectrum[0][centroid_low], mean=centroid_low, stddev = 75)
						fit_g_low = fitting.LevMarLSQFitter()
						g_low = fit_g_low(g_init_low, fit_channels_low, spectrum[0][fit_channels_low[0]:fit_channels_low[-1]+1])
						#print(fit_g.fit_info['param_cov'])
						#print(np.min(channelGrade[grade]))
						row[grade + 2] = (line_high-line_low)/(g_high.mean - g_low.mean)
						row[grade + 15] = (line_low*g_high.mean - line_high*g_low.mean)/(g_high.mean - g_low.mean)
						plt.figure()
						plt.plot(range(len(spectrum[0])), spectrum[0])
						plt.plot(fit_channels_high, g_high(fit_channels_high))
						plt.plot(fit_channels_low, g_low(fit_channels_low))
						plt.savefig('/Volumes/LaCie/CdTe/longGammaFlood/images/pixelGainOffset_Co/20170908_H100_long_gamma_Co57_-10_x' + str(x) + '_y' + str(y) + 'grade' + str(grade) + '_linefit.0V.eps')
						#plt.show()
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
t.writeto('/Volumes/LaCie/CdTe/longGammaFlood/20170908_H100_long_gamma_Co57_-10_gain_offset.fits')

