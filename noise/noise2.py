from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib as mpl
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

maxchannel = 1000
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
plt.imshow(countMap)
c = plt.colorbar()
c.set_label('Counts')
#plt.title(detector + ' ' + test + ' Pixel Map ' + '(' + etc + ')')
plt.tight_layout()
#plt.savefig('/disk/lif2/spike/detectorData/' + detector + '/figures/' + filename[:-4] + 'pixmap_corr.eps')
plt.savefig('/users/spike/det_figs/' + detector + '/' + filename[:-4] + 'pixmap_corr.eps')
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
			if g.fwhm < 1000:
				FWHM_map[row][col] = g.fwhm * gain[row][col]
			plt.hist(np.multiply(channelMap[row][col], gain[row][col]), bins = np.multiply(bins, gain[row][col]), range = (0-maxchannel*gain[row][col],maxchannel*gain[row][col]), histtype='stepfilled')
			plt.plot(np.multiply(fit_channels, gain[row][col]), g(fit_channels))
			plt.ylabel('Counts')
			if gainBool:
				plt.xlabel('Energy (keV)')
				plt.tight_layout()
				#plt.savefig('/disk/lif2/spike/detectorData/' + detector + '/figures/pixel_figs/' + filename[:-4] + 'x' + str(col) + 'y' + str(row) + '_spec_gain.eps')
				plt.savefig('/users/spike/det_figs/' + detector + '/pixels/' + filename[:-4] + 'x' + str(col) + 'y' + str(row) + '_spec_gain.eps')
			else:
				plt.xlabel('Channel')
				plt.tight_layout()
				#plt.savefig('/disk/lif2/spike/detectorData/' + detector + '/figures/pixel_figs/' + filename[:-4] + 'x' + str(col) + 'y' + str(row) + '_spec_corr.eps')
				plt.savefig('/users/spike/det_figs/' + detector + '/pixels/' + filename[:-4] + 'x' + str(col) + 'y' + str(row) + '_spec_corr.eps')
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
	plt.hist(FWHM, bins = 50, range = (0, 300), histtype='stepfilled')
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
	#plt.savefig('/disk/lif2/spike/detectorData/' + detector + '/figures/' + filename[:-4] + 'FWHMhist_gain.eps')
	plt.savefig('/users/spike/det_figs/' + detector + '/' + filename[:-4] + 'FWHMhist_gain.eps')
else:
	#plt.savefig('/disk/lif2/spike/detectorData/' + detector + '/figures/' + filename[:-4] + 'FWHMhist_corr.eps')
	plt.savefig('/users/spike/det_figs/' + detector + '/' + filename[:-4] + 'FWHMhist.eps')
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
masked = np.ma.masked_where(FWHM_map > 5, FWHM_map)
current_cmap = mpl.cm.get_cmap()
current_cmap.set_bad(color='gray')
plt.imshow(masked)
c = plt.colorbar()
if gainBool:
	c.set_label('FWHM (keV)')
else:
	c.set_label('FWHM (channels)')
#plt.title(detector + ' ' + test + ' FWHM Map ' + '(' + etc + ')')
plt.tight_layout()
if gainBool:
	#plt.savefig('/disk/lif2/spike/detectorData/' + detector + '/figures/' + filename[:-4] + 'FWHMmap_gain.eps')
	plt.savefig('/users/spike/det_figs/' + detector + '/' + filename[:-4] + 'FWHMmap_gain.eps')
else:
	#plt.savefig('/disk/lif2/spike/detectorData/' + detector + '/figures/' + filename[:-4] + 'FWHMmap_corr.eps')
	plt.savefig('/users/spike/det_figs/' + detector + '/' + filename[:-4] + 'FWHMmap.eps')
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
#plt.savefig('/disk/lif2/spike/detectorData/' + detector + '/figures/' + filename[:-4] + 'pixhist_corr.eps')
plt.savefig('/users/spike/det_figs/' + detector + '/' + filename[:-4] + 'pixhist.eps')
#plt.show()
plt.close()
'''
spectrum = np.histogram(data['PH'][START:END], bins = bins, range= (0-maxchannel, maxchannel))
plt.plot(spectrum[1][:-1], spectrum[0])
plt.xlabel('Channel')
plt.ylabel('Counts')
#plt.title(detector + ' ' + test + ' Spectrum ' + '(' + etc + ')')
plt.tight_layout()
plt.savefig('/disk/lif2/spike/detectorData/' + detector + '/figures/' + filename[:-4] + 'spec_corr.eps')
#plt.show()
plt.close()'''