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

if filepath[-1] == '/':
	filepath = filepath[:-1]

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

maxchannel = 250
bins = np.arange(0-maxchannel,maxchannel)

# There's no STIM data
channelMap = [[[[] for k in range(16)] for i in range(33)] for j in range(33)]
for i in np.arange(START, END):
	if data['UP'][i]:
		for j in range(9):
			channelMap[data['RAWY'][i] + (j//3) - 1][data['RAWX'][i] + (j%3) - 1][data['S_CAP'][i]].append(data['PH_RAW'][i][j])

file.close()

cap_offset = np.array([[[0 for c in range(16)] for x in range(32)] for y in range(32)])
for row in range(32):
	for col in range(32):
		for cap in range(16):

			plt.figure()
			#tempSpec = np.histogram(channelMap[row][col][cap], bins=bins, range = (0-maxchannel,maxchannel))
			tempSpec = plt.hist(channelMap[row][col][cap], bins=bins, range = (0-maxchannel,maxchannel))
			
			g_init = models.Gaussian1D(amplitude=np.max(tempSpec[0]), mean=0, stddev = 35)
			fit_g = fitting.LevMarLSQFitter()
			g = fit_g(g_init, tempSpec[1][:-1], tempSpec[0])
			cap_offset[row][col][cap] = g.mean.value
			plt.plot(tempSpec[1][:-1], g(tempSpec[1][:-1]))
			plt.text(g.mean.value + g.fwhm, g.amplitude.value, 'Mean: ' + str(round(g.mean.value, 2)))

			plt.xlabel('Channel')
			plt.ylabel('Counts')
			plt.tight_layout()
			plt.savefig('/users/spike/det_figs/' + detector + '/pixels/' + filename[:-5] + '_x' + str(col) + '_y' + str(row) + '_startcap_' + str(cap) + '.pdf')
			plt.close()

np.save(filepath[:-5] + '_startcap_offset.npy', cap_offset)