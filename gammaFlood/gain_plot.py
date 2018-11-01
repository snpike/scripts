import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os.path
import seaborn as sns
from matplotlib.gridspec import GridSpec
import scipy.stats as stats

sns.set_context('talk')
sns.set_style("ticks")
sns.set_palette("colorblind")

filepath = input('Please enter the filepath to the gain data: ').strip()

if filepath[-1] =='/':
	filepath = filepath[:-1]

detector = input('Please enter the detector ID: ').strip()

temp_response = input('Should I analyze the full detector? (y/n) ').strip()
while temp_response not in ['y','yes', 'n', 'no']:
    temp_response = input('Should I analyze the full detector? (y/n) ').strip()

full_det = (temp_response in ['y','yes'])
if full_det:
    region = [[0,31][0,31]]
else:
    region = [[int(input('Low column? ').strip()), int(input('Low row? ').strip())],[int(input('High column? ').strip()), int(input('High row? ').strip())]]

filename = filepath.split('/')[-1]

modifier = ''

filename_split = filename.split('.')

if len(filename_split) > 2:
	modifier = '.' + filename_split[-2]

gain = np.loadtxt(filepath)
gain_reg = gain[region[0][1]:region[1][1] + 1, region[0][0]:region[1][0] + 1]

plt.figure()
plt.imshow(gain_reg * 1000)
c = plt.colorbar()
c.set_label('Gain (eV/channel)')
plt.xticks(np.arange(0,region[1][0]-region[0][0]+1), [str(int(x)) for x in np.arange(region[0][0], region[1][0] + 1)])
plt.yticks(np.arange(0,region[1][1]-region[0][1]+1), [str(int(x)) for x in np.arange(region[0][1], region[1][1] + 1)])
plt.tight_layout()
plt.savefig('/users/spike/det_figs/' + detector + '/' + filename_split[0] + '_map' + modifier + '.pdf')
plt.close()