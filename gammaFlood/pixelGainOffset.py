from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from astropy.modeling import models, fitting

# From http://www.nndc.bnl.gov/nudat2/indx_dec.jsp
Co = 122.06
Am = 59.54
line_high = 10000
line_low = Am*line_high/Co

# Format of the input data
# events = np.array([[{'TIME': [], 'CHANNEL': [], 'GRADE': [], 'STIM': [], 'PH_COM': [], 'SOURCE': []} for i in range(32)] for j in range(32)])

# Columns: x, y, gain, offset
rows =[]

for x in range(32):
	for y in range(32):

		plt.figure()

		row = [x, y, 1.0, 0.0]

		dataCo = fits.open('/disk/lif2/spike/detectorData/H100/longGammaFlood/pixelData/H100_long_gamma_Co57_-10_0V_x' + str(x) + '_y' + str(y) + '.fits', memmap=True)[1].data
		channelCo = []
		for ch in dataCo['CHANNEL']:
			if ch < 20000:
				channelCo.append(ch)

		dataAm = fits.open('/disk/lif2/spike/detectorData/H100/longGammaFlood/pixelData/H100_long_gamma_Am241_-10_0V_x' + str(x) + '_y' + str(y) + '.fits', memmap=True)[1].data
		channelAm = []
		for ch in dataAm['CHANNEL']:
			if ch < 20000:
				channelAm.append(ch)

		if len(channelCo) and len(channelAm):
			spectrumCo = np.histogram(channelCo, bins=int(np.ceil(np.max(channelCo))))
			spectrumAm = np.histogram(channelAm, bins=int(np.ceil(np.max(channelAm))))

			plt.plot(range(len(spectrumCo[0])), spectrumCo[0], label = 'Co57')
			plt.plot(range(len(spectrumAm[0])), spectrumAm[0], label = 'Am241', color = 'r')
			plt.ylabel('Counts')
			plt.xlabel('Channel')
			plt.legend()

			if np.max(channelCo) > 12000 and np.max(channelAm) > 6000:
				centroid_high = np.argmax(spectrumCo[0][6000:12000]) + 6000
				fit_channels_high = np.arange(centroid_high-100, centroid_high + 250)
				g_init_high = models.Gaussian1D(amplitude=spectrumCo[0][centroid_high], mean=centroid_high, stddev = 75)
				fit_g_high = fitting.LevMarLSQFitter()
				g_high = fit_g_high(g_init_high, fit_channels_high, spectrumCo[0][fit_channels_high[0]:fit_channels_high[-1]+1])

				centroid_low = np.argmax(spectrumAm[0][2500:6000]) + 2500
				fit_channels_low = np.arange(centroid_low-100, centroid_low + 250)
				g_init_low = models.Gaussian1D(amplitude=spectrumAm[0][centroid_low], mean=centroid_low, stddev = 75)
				fit_g_low = fitting.LevMarLSQFitter()
				g_low = fit_g_low(g_init_low, fit_channels_low, spectrumAm[0][fit_channels_low[0]:fit_channels_low[-1]+1])
				
				row[2] = (line_high-line_low)/(g_high.mean - g_low.mean)
				row[3] = (line_low*g_high.mean - line_high*g_low.mean)/(g_high.mean - g_low.mean)
				plt.plot(fit_channels_high, g_high(fit_channels_high))
				plt.plot(fit_channels_low, g_low(fit_channels_low))
				
		plt.savefig('/disk/lif2/spike/detectorData/H100/figures/pixelFits/H100_long_gamma_Co57_Am241_-10_x' + str(x) + '_y' + str(y) + '_gain_offset_linefit.0V.eps')
		plt.close()

		rows.append(row)

columns = np.array(rows).T
fits_columns = []
col0 = fits.Column(name='RAWX', format='K', array=columns[0])
col1 = fits.Column(name='RAWY', format='K', array=columns[1])
fits_columns.append(col0)
fits_columns.append(col1)
fits_columns.append(fits.Column(name='GAIN', format='D', array=columns[2]))
fits_columns.append(fits.Column(name='OFFSET', format='D', array=columns[3]))
t = fits.BinTableHDU.from_columns(fits_columns)
t.writeto('/disk/lif2/spike/detectorData/H100/H100_long_gamma_Co57_Am241_-10_gain_offset.0V.fits')



