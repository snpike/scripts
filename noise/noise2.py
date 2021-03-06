from astropy.io import fits
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.modeling import models, fitting
import os.path
import pickle
import seaborn as sns
from matplotlib.gridspec import GridSpec

sns.set_context('talk')
sns.set_style("ticks")
sns.set_palette("colorblind")

filepath = input('Please enter the filepath to the noise data: ').strip()
while not os.path.exists(filepath):
	filepath = input('Please enter the filepath to the noise data: ').strip()

gainpath = input('Please enter the filepath to the gain data. If gain data is not available, press enter: ').strip()

slash = 0
i = 0
for char in filepath:
	if char == '/':
		slash = i
	i += 1

filename = filepath[slash + 1:]
detector = input('Please enter the detector ID: ').strip()
pos = int(input('What is the position of the detector? ').strip())
etc = input('Please enter any other important information (temperature, voltage, etc.): ')

gainBool = os.path.exists(gainpath)
gain = np.ones((32,32))
if gainBool:
	gain = pickle.load(open(gainpath, 'rb'))

test = 'noise'

file = fits.open(filepath)

data = file[1].data

mask = np.multiply((data['DET_ID'] == pos), (data['TEMP'] > -20))

START = np.argmax(mask)
END = len(mask) - np.argmax(mask[::-1])

maxchannel = 500
bins = np.arange(0-maxchannel,maxchannel)

# There's no STIM data
channelMap = [[[] for i in range(33)] for j in range(33)]
for i in np.arange(START, END):
	if data['UP'][i]:
		for j in range(9):
			channelMap[data['RAWY'][i] + (j//3) - 1][data['RAWX'][i] + (j%3) - 1].append(data['PH_RAW'][i][j])

'''		channelMap[data['RAWX'][i]-1][data['RAWY'][i]-1].append(data['PH_RAW'][i][0])
		channelMap[data['RAWX'][i]+0][data['RAWY'][i]-1].append(data['PH_RAW'][i][1])
		channelMap[data['RAWX'][i]+1][data['RAWY'][i]-1].append(data['PH_RAW'][i][2])
		channelMap[data['RAWX'][i]-1][data['RAWY'][i]+0].append(data['PH_RAW'][i][3])
		channelMap[data['RAWX'][i]+0][data['RAWY'][i]+0].append(data['PH_RAW'][i][4])
		channelMap[data['RAWX'][i]+1][data['RAWY'][i]+0].append(data['PH_RAW'][i][5])
		channelMap[data['RAWX'][i]-1][data['RAWY'][i]+1].append(data['PH_RAW'][i][6])
		channelMap[data['RAWX'][i]+0][data['RAWY'][i]+1].append(data['PH_RAW'][i][7])
		channelMap[data['RAWX'][i]+1][data['RAWY'][i]+1].append(data['PH_RAW'][i][8])'''

file.close()

countMap = [[len(channelMap[j][i]) for i in range(32)] for j in range(32)]
plt.figure()
plt.imshow(countMap, cmap='coolwarm')
c = plt.colorbar()
c.set_label('Counts')
#plt.title(detector + ' ' + test + ' Pixel Map ' + '(' + etc + ')')
plt.tight_layout()
#plt.savefig('/disk/lif2/spike/detectorData/' + detector + '/figures/' + filename[:-4] + 'pixmap_corr.pdf')
plt.savefig('/disk/lif2/spike/det_figs/' + detector + '/' + filename[:-4] + 'pixmap.pdf')
#plt.show()
plt.close()

FWHM = []
FWHM_map = np.array([[np.nan for i in range(32)] for j in range(32)])
for row in range(32):
	for col in range(32):
		if (channelMap[row][col]):
			tempSpec = np.histogram(channelMap[row][col], bins=bins, range = (0-maxchannel,maxchannel))
			#centroid = np.argmax(tempSpec[0])
			fit_channels = tempSpec[1][:-1]
			g_init = models.Gaussian1D(amplitude=np.max(tempSpec[0]), mean=0, stddev = 75)
			fit_g = fitting.LevMarLSQFitter()
			g = fit_g(g_init, fit_channels, tempSpec[0])
			FWHM.append(g.fwhm * gain[row][col])
			FWHM_map[row][col] = g.fwhm * gain[row][col]
			plt.hist(np.multiply(channelMap[row][col], gain[row][col]), bins = np.multiply(bins, gain[row][col]), range = (0-maxchannel*gain[row][col],maxchannel*gain[row][col]), histtype='stepfilled')
			plt.plot(np.multiply(fit_channels, gain[row][col]), g(fit_channels))
			plt.ylabel('Counts')
			if gainBool:
				plt.xlabel('Energy (keV)')
				plt.text(gain[row][col] * (g.mean + g.fwhm/2), g.amplitude, 'Mean = ' + str(gain[row][col] * g.mean) + ' keV', fontsize = 12)
				plt.text(gain[row][col] * (g.mean + g.fwhm/2), g.amplitude*0.8, 'FWHM = ' + str(gain[row][col] * g.fwhm) + ' keV', fontsize = 12)
				plt.tight_layout()
				#plt.savefig('/disk/lif2/spike/detectorData/' + detector + '/figures/pixel_figs/' + filename[:-4] + 'x' + str(col) + 'y' + str(row) + '_spec_gain.pdf')
				plt.savefig('/disk/lif2/spike/det_figs/' + detector + '/pixels/' + filename[:-4] + 'x' + str(col) + 'y' + str(row) + '_spec_gain.pdf')
			else:
				plt.xlabel('Channel')
				plt.text(g.mean.value + g.fwhm/2, g.amplitude.value, 'Mean = ' + str(int(round(g.mean.value, 0))) + ' channels', fontsize = 12)
				plt.text(g.mean.value + g.fwhm/2, g.amplitude.value*0.8, 'FWHM = ' + str(int(round(g.fwhm, 0))) + ' channels', fontsize = 12)
				plt.tight_layout()
				#plt.savefig('/disk/lif2/spike/detectorData/' + detector + '/figures/pixel_figs/' + filename[:-4] + 'x' + str(col) + 'y' + str(row) + '_spec_corr.pdf')
				plt.savefig('/disk/lif2/spike/det_figs/' + detector + '/pixels/' + filename[:-4] + 'x' + str(col) + 'y' + str(row) + '_spec.pdf')
			plt.close()

#FWHM_hist = np.histogram(FWHM, bins = 50, range = (0, 300))
plt.figure()
if gainBool:
	plt.hist(FWHM, bins = 50, range = (0, 4), histtype='stepfilled')
	bot, top = plt.ylim()
	left, right = plt.xlim()
	plt.text(right*0.5, top*0.8, 'Mean = ' + str(int(round(np.mean(FWHM) * 1000, 0))) + ' eV', fontsize = 16)
	plt.text(right*0.5, top*0.6, '1-Sigma = ' + str(int(round(np.std(FWHM) * 1000, 0))) + ' eV', fontsize = 16)
	plt.xlabel('FWHM (keV)')
else:
	plt.hist(FWHM, bins = 50, range = (0, 150), histtype='stepfilled')
	bot, top = plt.ylim()
	left, right = plt.xlim()
	plt.text(right*0.5, top*0.8, 'Mean = ' + str(round(np.mean(FWHM), 0)) + ' channels', fontsize = 16)
	plt.text(right*0.5, top*0.6, '1-Sigma = ' + str(round(np.std(FWHM), 0)) + ' channels', fontsize = 16)
	plt.xlabel('FWHM (channels)')
plt.ylabel('Pixels')
#plt.xlim((0, 300))
#plt.ylim(ymin=0)
#plt.title(detector + ' ' + test + ' FWHM Histogram ' + '(' + etc + ')')
plt.tight_layout()
if gainBool:
	#plt.savefig('/disk/lif2/spike/detectorData/' + detector + '/figures/' + filename[:-4] + 'FWHMhist_gain.pdf')
	plt.savefig('/disk/lif2/spike/det_figs/' + detector + '/' + filename[:-4] + 'FWHMhist_gain.pdf')
else:
	#plt.savefig('/disk/lif2/spike/detectorData/' + detector + '/figures/' + filename[:-4] + 'FWHMhist_corr.pdf')
	plt.savefig('/disk/lif2/spike/det_figs/' + detector + '/' + filename[:-4] + 'FWHMhist.pdf')
#plt.show()
plt.close()

outfile = open('/disk/lif2/spike/detectorData/' + detector + '/noise2.out', 'w')
outfile.write('Mean FWHM: '  +'\n')
outfile.write(str(np.mean(FWHM)) + '\n')
outfile.write('Std dev FWHM: '  +'\n')
outfile.write(str(np.std(FWHM)) + '\n')
outfile.close()

dumpfile = open('/disk/lif2/spike/detectorData/' + detector + '/noise2_FWHMmap_dump.txt', 'wb')
pickle.dump(FWHM_map, dumpfile)
dumpfile.close()

plt.figure()
if gainBool:
	masked = np.ma.masked_where(FWHM_map > 5, FWHM_map)
else:
	masked = np.ma.masked_where(FWHM_map > 500, FWHM_map)
current_cmap = mpl.cm.get_cmap()
current_cmap.set_bad(color='gray')
plt.imshow(masked, cmap='coolwarm')
c = plt.colorbar()
if gainBool:
	c.set_label('FWHM (keV)')
else:
	c.set_label('FWHM (channels)')
#plt.title(detector + ' ' + test + ' FWHM Map ' + '(' + etc + ')')
plt.tight_layout()
if gainBool:
	#plt.savefig('/disk/lif2/spike/detectorData/' + detector + '/figures/' + filename[:-4] + 'FWHMmap_gain.pdf')
	plt.savefig('/disk/lif2/spike/det_figs/' + detector + '/' + filename[:-4] + 'FWHMmap_gain.pdf')
else:
	#plt.savefig('/disk/lif2/spike/detectorData/' + detector + '/figures/' + filename[:-4] + 'FWHMmap_corr.pdf')
	plt.savefig('/disk/lif2/spike/det_figs/' + detector + '/' + filename[:-4] + 'FWHMmap.pdf')
#plt.show()
plt.close()


#noiseHist = np.histogram(np.array(countMap).flatten(), bins = 50)
plt.figure()
plt.hist(np.array(countMap).flatten(), bins = 50, histtype = 'stepfilled')
# plt.xlim(noiseHist[1][0], noiseHist[1][-1])
# plt.ylim(ymin=0)
plt.ylabel('Pixels')
plt.xlabel('Counts')
#plt.xticks(noiseHist[1][1:-1])
#plt.title(detector + ' ' + test + ' Count Histogram ' + '(' + etc + ')')
plt.tight_layout()
#plt.savefig('/disk/lif2/spike/detectorData/' + detector + '/figures/' + filename[:-4] + 'pixhist_corr.pdf')
plt.savefig('/disk/lif2/spike/det_figs/' + detector + '/' + filename[:-4] + 'pixhist.pdf')
#plt.show()
plt.close()
'''
spectrum = np.histogram(data['PH'][START:END], bins = bins, range= (0-maxchannel, maxchannel))
plt.plot(spectrum[1][:-1], spectrum[0])
plt.xlabel('Channel')
plt.ylabel('Counts')
#plt.title(detector + ' ' + test + ' Spectrum ' + '(' + etc + ')')
plt.tight_layout()
plt.savefig('/disk/lif2/spike/detectorData/' + detector + '/figures/' + filename[:-4] + 'spec_corr.pdf')
#plt.show()
plt.close()'''