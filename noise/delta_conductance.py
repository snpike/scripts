import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_context('talk')
sns.set_style('ticks')
sns.set_palette('colorblind')

filepath_i = input('Please enter the path to the initial interpixel data: ').strip()
filepath_f = input('Please enter the path to the initial interpixel data: ').strip()
date_i = (filepath_i.split('/')[-3]).split('_')[0]
date_f = (filepath_f.split('/')[-3]).split('_')[0]
detector = input('Please enter the detector ID: ').strip()

row_i = []
col_i = []
row_f = []
col_f = []

i=0

file = open(filepath_i, 'r')
for line in file:
	if 4<i<69:
		if not i%2:
			row_i.append(line.split()[1:])
	if 71<i<137:
		if not i%2:
			col_i.append(line.split()[1:])
	i+=1

i=0

file = open(filepath_f, 'r')
for line in file:
	if 4<i<69:
		if not i%2:
			row_f.append(line.split()[1:])
	if 71<i<137:
		if not i%2:
			col_f.append(line.split()[1:])
	i+=1

row_i = np.array(row_i).astype(float)
row_f = np.array(row_f).astype(float)
col_i = np.array(col_i).astype(float).T
col_f = np.array(col_f).astype(float).T

row_diff = row_f - row_i
col_diff = col_f - col_i

vmax = np.max([np.max(row_diff), np.max(col_diff)])
plt.figure()
plt.imshow(row_diff, vmin=0, vmax=vmax, cmap='coolwarm')
c = plt.colorbar()
c.set_label('(S/1.7e10)')
plt.tight_layout()
plt.savefig('/users/spike/det_figs/' + detector + '/' + date_f +'_m_' + date_i + '_deltaconductance' + '_row_map.pdf')
plt.close()

plt.figure()
plt.imshow(col_diff, vmin=0, vmax=vmax, cmap='coolwarm')
c = plt.colorbar()
c.set_label('(S/1.7e10)')
plt.tight_layout()
plt.savefig('/users/spike/det_figs/' + detector + '/' + date_f +'_m_' + date_i + '_deltaconductance' + '_col_map.pdf')
plt.close()