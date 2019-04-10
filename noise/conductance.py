import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_context('talk')
sns.set_style('ticks')
sns.set_palette('RdBu_r')

filepath = input('Please enter the path to the interpixel data: ').strip()
pathsplit = filepath.split('/')
detector = input('Please enter the detector ID: ').strip()

row_data = []
col_data = []

i=0

file = open(filepath, 'r')
for line in file:
	if 4<i<69:
		if not i%2:
			row_data.append(line.split()[1:])
	if 71<i<137:
		if not i%2:
			col_data.append(line.split()[1:])
	i+=1

row_data = np.array(row_data).astype(int)
col_data = np.array(col_data).astype(int).T

vmax = np.max([np.max(row_data), np.max(col_data)])
plt.figure()
plt.imshow(row_data, vmin=0, vmax=vmax)
c = plt.colorbar()
c.set_label('(S/1.7e10)')
plt.tight_layout()
plt.savefig('/users/spike/det_figs/' + detector + '/' + pathsplit[-3] + '_row_map.pdf')
plt.close()

plt.figure()
plt.imshow(col_data, vmin=0, vmax=vmax)
c = plt.colorbar()
c.set_label('(S/1.7e10)')
plt.tight_layout()
plt.savefig('/users/spike/det_figs/' + detector + '/' + pathsplit[-3] + '_col_map.pdf')
plt.close()