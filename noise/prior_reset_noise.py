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
priorMap = [[[] for i in range(33)] for j in range(33)]
resetMap = [[[] for i in range(33)] for j in range(33)]
time_riseMap = [[[] for i in range(33)] for j in range(33)]
channel = []
prior = []
reset = []
time_rise = []
for i in np.arange(START, END):
	if data['UP'][i]:
		for j in range(9):
			channelMap[data['RAWY'][i] + (j//3) - 1][data['RAWX'][i] + (j%3) - 1].append(data['PH_RAW'][i][j])
			priorMap[data['RAWY'][i] + (j//3) - 1][data['RAWX'][i] + (j%3) - 1].append(data['PRIOR'][i])
			resetMap[data['RAWY'][i] + (j//3) - 1][data['RAWX'][i] + (j%3) - 1].append(data['RESET'][i])
			time_riseMap[data['RAWY'][i] + (j//3) - 1][data['RAWX'][i] + (j%3) - 1].append(data['NUMRISE'][i]/data['DENRISE'][i])
			channel.append(data['PH_RAW'][i][j])
			prior.append(data['PRIOR'][i])
			reset.append(data['RESET'][i])
			time_rise.append(data['NUMRISE'][i]/data['DENRISE'][i])


prior = np.array(prior)
channel = np.array(channel)
reset = np.array(reset)
rand = np.random.randint(0, len(prior), size = 10000)
plt.figure()
plt.scatter(prior[rand], channel[rand], s=1, marker = '.', rasterized=True)
plt.xlabel('Time since last event')
plt.ylabel('Channel')
plt.tight_layout()
plt.savefig('/users/spike/det_figs/' + detector + '/' + filename[:-4] + 'last_event.pdf')
plt.close()

plt.figure()
plt.hist(prior[rand], bins = 50)
plt.xlabel('Time since last event')
plt.ylabel('Counts')
plt.tight_layout()
plt.savefig('/users/spike/det_figs/' + detector + '/' + filename[:-4] + 'last_event_hist.pdf')
plt.close()

plt.figure()
plt.scatter(reset[rand], channel[rand], s = 1, marker = '.', rasterized=True)
plt.xlabel('Time since last reset')
plt.ylabel('Channel')
plt.tight_layout()
plt.savefig('/users/spike/det_figs/' + detector + '/' + filename[:-4] + 'last_reset.pdf')
plt.close()

plt.figure()
plt.hist(reset[rand], bins = 50)
plt.xlabel('Time since last reset')
plt.ylabel('Counts')
plt.tight_layout()
plt.savefig('/users/spike/det_figs/' + detector + '/' + filename[:-4] + 'last_reset_hist.pdf')
plt.close()

plt.figure()
plt.scatter(time_rise[rand], channel[rand], s = 1, marker = '.', rasterized=True)
plt.xlabel('Time of rise conversion')
plt.ylabel('Channel')
plt.tight_layout()
plt.savefig('/users/spike/det_figs/' + detector + '/' + filename[:-4] + 'time_rise.pdf')
plt.close()

plt.figure()
plt.hist(time_rise[rand], bins = 50)
plt.xlabel('Time of rise conversion')
plt.ylabel('Counts')
plt.tight_layout()
plt.savefig('/users/spike/det_figs/' + detector + '/' + filename[:-4] + 'time_rise_hist.pdf')
plt.close()

for row in range(32):
	for col in range(32):
		plt.figure()
		plt.scatter(priorMap[col][row], channelMap[col][row], s=1, marker = '.', rasterized=True)
		plt.xlabel('Time since last event')
		plt.ylabel('Channel')
		plt.tight_layout()
		plt.savefig('/users/spike/det_figs/' + detector + '/pixels/' + filename[:-4] + 'last_event_x' + str(col) + '_y' + str(row) + '.pdf')
		plt.close()

		plt.figure()
		plt.hist(priorMap[col][row], bins = 50)
		plt.xlabel('Time since last event')
		plt.ylabel('Counts')
		plt.tight_layout()
		plt.savefig('/users/spike/det_figs/' + detector + '/pixels/' + filename[:-4] + 'last_event_hist_x' + str(col) + '_y' + str(row) + '.pdf')
		plt.close()

		plt.figure()
		plt.scatter(resetMap[col][row], channelMap[col][row], s = 1, marker = '.', rasterized=True)
		plt.xlabel('Time since last reset')
		plt.ylabel('Channel')
		plt.tight_layout()
		plt.savefig('/users/spike/det_figs/' + detector + '/pixels/' + filename[:-4] + 'last_reset_x' + str(col) + '_y' + str(row) + '.pdf')
		plt.close()

		plt.figure()
		plt.hist(resetMap[col][row], bins = 50)
		plt.xlabel('Time since last reset')
		plt.ylabel('Counts')
		plt.tight_layout()
		plt.savefig('/users/spike/det_figs/' + detector + '/pixels/' + filename[:-4] + 'last_reset_hist_x' + str(col) + '_y' + str(row) + '.pdf')
		plt.close()

		plt.figure()
		plt.scatter(time_riseMap[col][row], channelMap[col][row], s = 1, marker = '.', rasterized=True)
		plt.xlabel('Time of rise conversion')
		plt.ylabel('Channel')
		plt.tight_layout()
		plt.savefig('/users/spike/det_figs/' + detector + '/pixels/' + filename[:-4] + 'time_rise_x' + str(col) + '_y' + str(row) + '.pdf')
		plt.close()

		plt.figure()
		plt.hist(time_riseMap[col][row], bins = 50)
		plt.xlabel('Time of rise conversion')
		plt.ylabel('Counts')
		plt.tight_layout()
		plt.savefig('/users/spike/det_figs/' + detector + '/pixels/' + filename[:-4] + 'time_rise_hist_x' + str(col) + '_y' + str(row) + '.pdf')
		plt.close()



