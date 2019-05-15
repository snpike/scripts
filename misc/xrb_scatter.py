import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.gridspec import GridSpec
from scipy.optimize import curve_fit
import csv

sns.set_context('talk', font_scale = 1.5)
sns.set_style("ticks")
sns.set_palette("colorblind")


plt.figure(figsize=(12,8))

marker_dict = {'LMXB': 'ko', 'SGXB': 'ks', 'Be transient': 'kx', 'ULXP': 'c*'}
markersize_dict = {'LMXB': 10, 'SGXB': 10, 'Be transient': 10, 'ULXP': 25}

with open('/Users/sean/Desktop/Bildsten XRBs.csv', 'r') as file:
	reader = csv.DictReader(file)
	for row in reader:
		if row['Porb'] != '999':
			if float(row['Porb']) > 1.0:
				plt.plot(float(row['Porb']), float(row['Pspin']), marker_dict[row['Type']], markersize=markersize_dict[row['Type']])
		# else:
		# 	plt.axhline(float(row['Pspin']), color='cyan', linewidth = 5)

plt.loglog()
plt.ylabel('Spin Period (s)')
plt.xlabel('Orbital Period (days)')
plt.tight_layout()
plt.savefig('/Users/sean/Desktop/bildsten.pdf')
plt.close()

