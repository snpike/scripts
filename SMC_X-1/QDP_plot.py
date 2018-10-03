import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_context('paper', font_scale=1.2, rc={"lines.linewidth": 0.8})
sns.set_style("ticks")
sns.set_palette("colorblind")

'''
pl_ratio = ['/Users/sean/Desktop/SMC_X-1_final/products001_part1/pl_ratio_epoch1.qdp', \
			'/Users/sean/Desktop/SMC_X-1_final/products001_part2/pl_ratio_epoch2.qdp', \
			'/Users/sean/Desktop/SMC_X-1_final/products003/pl_ratio_epoch3.qdp']

cutoffpl_spectra = ['/Users/sean/Desktop/SMC_X-1_final/products001_part1/tbabs_tbpcf_cutoffpl_gauss_gauss.qdp', \
					'/Users/sean/Desktop/SMC_X-1_final/products001_part2/tbabs_tbpcf_cutoffpl_gauss_gauss.qdp', \
					'/Users/sean/Desktop/SMC_X-1_final/products003/tbabs_tbpcf_cutoffpl_gauss_gauss.qdp']

cutoffpl_ratio = ['/Users/sean/Desktop/SMC_X-1_final/products001_part1/tbabs_tbpcf_cutoffpl_gauss_gauss_ratio.qdp', \
				  '/Users/sean/Desktop/SMC_X-1_final/products001_part2/tbabs_tbpcf_cutoffpl_gauss_gauss_ratio.qdp', \
				  '/Users/sean/Desktop/SMC_X-1_final/products003/tbabs_tbpcf_cutoffpl_gauss_gauss_ratio.qdp']
'''

epochs = {'epochI': ['/Users/sean/Desktop/SMC_X-1_final/products001_part1/eeufspec_const_epoch1.qdp', \
		  '/Users/sean/Desktop/SMC_X-1_final/products001_part1/pl_ratio_epoch1.qdp', \
		  '/Users/sean/Desktop/SMC_X-1_final/products001_part1/tbpcf_cutoffpl_gauss_gauss_frozengamma_ratio.qdp', \
		  '/Users/sean/Desktop/SMC_X-1_final/products001_part1/tbabs_fdcut_bb_gauss_gauss.qdp'], \
		  'epochII': ['/Users/sean/Desktop/SMC_X-1_final/products001_part2/eeufspec_const_epoch2.qdp', \
		   '/Users/sean/Desktop/SMC_X-1_final/products001_part2/pl_ratio_epoch2.qdp', \
		   '/Users/sean/Desktop/SMC_X-1_final/products001_part2/tbpcf_cutoffpl_gauss_gauss_frozengamma_ratio.qdp', \
		   '/Users/sean/Desktop/SMC_X-1_final/products001_part2/tbabs_fdcut_bb_gauss_gauss_ratio.qdp'], \
		   'epochIII': ['/Users/sean/Desktop/SMC_X-1_final/products003/eeufspec_const_epoch3.qdp', \
			'/Users/sean/Desktop/SMC_X-1_final/products003/pl_ratio_epoch3.qdp', \
			'/Users/sean/Desktop/SMC_X-1_final/products003/tbpcf_cutoffpl_bb_bb_gauss_frozengamma_ratio.qdp', \
			'/Users/sean/Desktop/SMC_X-1_final/products003/fdcut_bb_bb_gauss_ratio.qdp']}

for key in epochs:
	file = open(epochs[key][0], 'r')

	data = [[],[]]
	i = 0
	j = 0
	for line in file:
		if i >2:
			temp = line.split()
			if temp[0] != 'NO':
				data[j].append(temp)
			else:
				j += 1
		i += 1

	data[0] = np.array(data[0]).astype(float).T
	data[1] = np.array(data[1]).astype(float).T

	fig, (ax1, ax2, ax3, ax4) = plt.subplots(4,1, sharex=True,  gridspec_kw = {'height_ratios':[3, 1, 1, 1], 'hspace':0})
	ax1.errorbar(data[0][0],data[0][2], xerr=data[0][1], yerr=data[0][3], ls = 'none', color='black', label='FPMA')	#spectrum
	if key=='epochIII':
		ax1.set_ylim((0.205, np.max(data[0][2]) * 1.2))
	else:
		ax1.set_ylim((0.0105, 0.1))
	# ax1.plot(data[0][0],data[0][4], ls = 'steps-mid', color = 'black')	#total model?
	# ax1.plot(data[0][0],data[0][5], ls = '--', color = 'black')	#cutoffpl
	# ax1.plot(data[0][0],data[0][6], ls = '--', color = 'black')	#13keV bump
	# ax1.plot(data[0][0],data[0][7], ls = '--', color = 'black')	#Fe line
	ax1.errorbar(data[1][0],data[1][2], xerr=data[1][1], yerr=data[1][3], ls = 'none', color='red', label='FPMB')	#spectrum
	# ax1.plot(data[1][0],data[1][4], ls = 'steps-mid', color = 'red')	#total model?
	# ax1.plot(data[1][0],data[1][5], ls = '--', color = 'red')	#cutoffpl
	# ax1.plot(data[1][0],data[1][6], ls = '--', color = 'red')	#13keV bump
	# ax1.plot(data[1][0],data[1][7], ls = '--', color = 'red')	#Fe line
	ax1.set_xscale('log')
	ax1.set_yscale('log')
	ax1.set_xlim((3,40))
	ax1.set_ylabel(r'$\mathrm{keV^2(\gamma\ cm^{-2}\ s^{-1}\ keV^{-1})}$')
	ax1.legend(loc=1)
	#plt.show()
	#plt.close()

	file = open(epochs[key][1], 'r')

	data = [[],[]]
	i = 0
	j = 0
	for line in file:
		if i >2:
			temp = line.split()
			if temp[0] != 'NO':
				data[j].append(temp)
			else:
				j += 1
		i += 1

	data[0] = np.array(data[0]).astype(float).T
	data[1] = np.array(data[1]).astype(float).T

	ax2.errorbar(data[0][0],data[0][2]-1, xerr=data[0][1], yerr=data[0][3], ls = 'none', color='black')	#pl ratio
	ax2.errorbar(data[1][0],data[1][2]-1, xerr=data[1][1], yerr=data[1][3], ls = 'none', color='red')
	ax2.text(3.1, -0.7, r'$tbabs \times powerlaw$')
	ax2.axhline(y=0, lw = 0.8)
	ax2.set_ylim((-0.8,0.8))
	ax2.set_yscale('linear')

	file = open(epochs[key][3], 'r')

	data = [[],[]]
	i = 0
	j = 0
	for line in file:
		if i >2:
			temp = line.split()
			if temp[0] != 'NO':
				data[j].append(temp)
			else:
				j += 1
		i += 1

	data[0] = np.array(data[0]).astype(float).T
	data[1] = np.array(data[1]).astype(float).T

	ax3.errorbar(data[0][0],data[0][2]-1, xerr=data[0][1], yerr=data[0][3], ls = 'none', color='black')	#fdcut ratio
	ax3.errorbar(data[1][0],data[1][2]-1, xerr=data[1][1], yerr=data[1][3], ls = 'none', color='red')
	ax3.axhline(y=0, lw = 0.8)
	ax3.set_ylim((-0.15, 0.15))
	if key == 'epochIII':
		ax3.text(3.1, -0.13125, r'$tbabs \times (fdcut + bbody_{0.3} + bbody_{1.5} + gauss_{6.4})$')
	else:
		ax3.text(3.1, -0.13125, r'$tbabs \times (fdcut + bbody_{0.3} + gauss_{6.4} + gauss_{13.5})$')
	ax3.set_yscale('linear')
	ax3.set_ylabel('(data-model)/model')

	file = open(epochs[key][2], 'r')

	data = [[],[]]
	i = 0
	j = 0
	for line in file:
		if i >2:
			temp = line.split()
			if temp[0] != 'NO':
				data[j].append(temp)
			else:
				j += 1
		i += 1

	data[0] = np.array(data[0]).astype(float).T
	data[1] = np.array(data[1]).astype(float).T

	ax4.errorbar(data[0][0],data[0][2]-1, xerr=data[0][1], yerr=data[0][3], ls = 'none', color='black')	#cutoffpl ratio
	ax4.errorbar(data[1][0],data[1][2]-1, xerr=data[1][1], yerr=data[1][3], ls = 'none', color='red')
	ax4.axhline(y=0, lw = 0.8)
	ax4.set_ylim((-0.15, 0.15))
	ax4.set_yscale('linear')
	ax4.set_xlabel('Energy (keV)')
	plt.xticks([3,4,5,6,7,8,9,10,20,30, 40], [3,4,5,6,7,8,9,10,20,30, 40])
	ax4.set_ylim((-0.15,0.15))
	if key == 'epochIII':
		ax4.text(3.1, -0.13125, r'$tbpcf \times (cutoffpl + bbody_{0.3} + bbody_{1.5} + gauss_{6.4})$')
	else:
		ax4.text(3.1, -0.13125, r'$tbpcf \times (cutoffpl + gauss_{6.4} + gauss_{13.5})$')
	fig.tight_layout(h_pad = 0)
	#plt.show()
	plt.savefig('/Users/sean/Desktop/MyPapers/SMC_X-1_letter/figures/spectral_analysis_' + key + '.pdf')
	plt.close()