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
channel = []
prior = []
reset = []
for i in np.arange(START, END):
	if data['UP'][i]:
		for j in range(9):
			channelMap[data['RAWY'][i] + (j//3) - 1][data['RAWX'][i] + (j%3) - 1].append(data['PH_RAW'][i][j])
			channel.append(data['PH_RAW'][i][j])
			prior.append(data['PRIOR'][i])
			reset.append(data['RESET'][i])

plt.figure()
plt.scatter(prior, channel, marker = '.')
plt.show()
plt.close()

plt.figure()
plt.scatter(reset, channel, marker = '.')
plt.show()
plt.close()

