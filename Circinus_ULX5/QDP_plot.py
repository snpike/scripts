import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_context('poster')
sns.set_style("ticks")
sns.set_palette("colorblind")




file = open('/Users/sean/Desktop/HEAD2019/tbabs_diskbb_frozennH_xrt_nustar_10keV_eeuf.qdp', 'r')

data = [[] for x in range(6)]
i = 0
j = 0
for line in file:
	if i >2:
		temp = line.split()
		if temp[0] != 'NO':
			if temp[-1]=='NO':
				temp = temp[:-1]
			data[j].append(temp)
		else:
			j += 1
	i += 1

labels = ['Swift XRT', 'NuSTAR FPMA', 'NuSTAR FPMB']

for i in range(6):
	data[i] = np.array(data[i]).astype(float).T

fig, (ax1, ax2) = plt.subplots(2,1, sharex=True,  gridspec_kw = {'height_ratios':[3, 1], 'hspace':0}, figsize=(12, 8))
for i in range(3):
	ax1.errorbar(data[i][0],data[i][2], xerr=data[i][1], yerr=data[i][3], ls = 'none', color='C' + str(i), label=labels[i])
	ax1.plot(data[i][0],data[i][4], ls = 'steps-mid', color = 'C' + str(i))	#total model
	ax1.set_xscale('log')
	ax1.set_yscale('log')
	# ax1.set_xlim((3,40))
ax1.set_ylabel(r'$\mathrm{keV^2(\gamma\ cm^{-2}\ s^{-1}\ keV^{-1})}$')
ax1.legend()

for i in range(3):
	j = i+3
	ax2.errorbar(data[j][0],data[j][2], xerr=data[j][1], yerr=data[j][3], ls = 'none', color='C' + str(i))	#pl ratio
	ax2.axhline(y=1, lw = 0.8)
	# ax2.set_ylim((-0.8,0.8))
	ax2.set_yscale('log')
ax2.set_ylabel('data/model')
ax2.set_xlabel('Energy (keV)')

plt.tight_layout()
plt.savefig('/Users/sean/Desktop/HEAD2019/tbabs_diskbb_frozennH_xrt_nustar_10keV_eeuf.pdf')
plt.close()


file = open('/Users/sean/Desktop/HEAD2019/tbabs_diskbb_frozennH_xrt_nustar_30keV_eeuf.qdp', 'r')

data = [[] for x in range(6)]
i = 0
j = 0
for line in file:
	if i >2:
		temp = line.split()
		if temp[0] != 'NO':
			if temp[-1]=='NO':
				temp = temp[:-1]
			data[j].append(temp)
		else:
			j += 1
	i += 1

labels = ['Swift XRT', 'NuSTAR FPMA', 'NuSTAR FPMB']

for i in range(6):
	data[i] = np.array(data[i]).astype(float).T

fig, (ax1, ax2) = plt.subplots(2,1, sharex=True,  gridspec_kw = {'height_ratios':[3, 1], 'hspace':0}, figsize=(12, 8))
for i in range(3):
	ax1.errorbar(data[i][0],data[i][2], xerr=data[i][1], yerr=data[i][3], ls = 'none', color='C' + str(i), label=labels[i])
	ax1.plot(data[i][0],data[i][4], ls = 'steps-mid', color = 'C' + str(i))	#total model
	ax1.set_xscale('log')
	ax1.set_yscale('log')
	# ax1.set_xlim((3,40))
ax1.set_ylabel(r'$\mathrm{keV^2(\gamma\ cm^{-2}\ s^{-1}\ keV^{-1})}$')
ax1.legend()

for i in range(3):
	j = i+3
	ax2.errorbar(data[j][0],data[j][2], xerr=data[j][1], yerr=data[j][3], ls = 'none', color='C' + str(i))	#pl ratio
	ax2.axhline(y=1, lw = 0.8)
	# ax2.set_ylim((-0.8,0.8))
	ax2.set_yscale('log')
ax2.set_ylabel('data/model')
ax2.set_yticks([1, 100])
ax2.set_xlabel('Energy (keV)')

plt.tight_layout()
plt.savefig('/Users/sean/Desktop/HEAD2019/tbabs_diskbb_frozennH_xrt_nustar_30keV_eeuf.pdf')
plt.close()

file = open('/Users/sean/Desktop/HEAD2019/tbabs_diskbb_frozennH_xrt_nustar_10keV_ldata.qdp', 'r')

data = [[] for x in range(6)]
i = 0
j = 0
for line in file:
	if i >2:
		temp = line.split()
		if temp[0] != 'NO':
			while temp[-1]=='NO':
				temp = temp[:-1]
			data[j].append(temp)
		else:
			j += 1
	i += 1

labels = ['Swift XRT', 'NuSTAR FPMA', 'NuSTAR FPMB']

for i in range(6):
	data[i] = np.array(data[i]).astype(float).T

fig, (ax1, ax2) = plt.subplots(2,1, sharex=True,  gridspec_kw = {'height_ratios':[3, 1], 'hspace':0}, figsize=(12, 8))
for i in range(3):
	ax1.errorbar(data[i][0],data[i][2], xerr=data[i][1], yerr=data[i][3], ls = 'none', color='C' + str(i), label=labels[i])
	ax1.plot(data[i][0],data[i][-1], ls = 'steps-mid', color = 'C' + str(i))	#total model
	ax1.set_xscale('log')
	ax1.set_yscale('log')
	# ax1.set_xlim((3,40))
ax1.set_ylabel(r'$\mathrm{\gamma\ cm^{-2}\ s^{-1}\ keV^{-1}}$')
ax1.legend()

for i in range(3):
	j = i+3
	ax2.errorbar(data[j][0],data[j][2], xerr=data[j][1], yerr=data[j][3], ls = 'none', color='C' + str(i))	#pl ratio
	ax2.axhline(y=1, lw = 0.8)
	# ax2.set_ylim((-0.8,0.8))
	ax2.set_yscale('log')
ax2.set_ylabel('data/model')
ax2.set_xlabel('Energy (keV)')

plt.tight_layout()
plt.savefig('/Users/sean/Desktop/HEAD2019/tbabs_diskbb_frozennH_xrt_nustar_10keV_ldata.pdf')
plt.close()

file = open('/Users/sean/Desktop/HEAD2019/tbabs_diskbb_frozennH_xrt_nustar_30keV_ldata.qdp', 'r')

data = [[] for x in range(6)]
i = 0
j = 0
for line in file:
	if i >2:
		temp = line.split()
		if temp[0] != 'NO':
			while temp[-1]=='NO':
				temp = temp[:-1]
			data[j].append(temp)
		else:
			j += 1
	i += 1

labels = ['Swift XRT', 'NuSTAR FPMA', 'NuSTAR FPMB']

for i in range(6):
	data[i] = np.array(data[i]).astype(float).T

fig, (ax1, ax2) = plt.subplots(2,1, sharex=True,  gridspec_kw = {'height_ratios':[3, 1], 'hspace':0}, figsize=(12, 8))
for i in range(3):
	ax1.errorbar(data[i][0],data[i][2], xerr=data[i][1], yerr=data[i][3], ls = 'none', color='C' + str(i), label=labels[i])
	ax1.plot(data[i][0],data[i][-1], ls = 'steps-mid', color = 'C' + str(i))	#total model
	ax1.set_xscale('log')
	ax1.set_yscale('log')
	# ax1.set_xlim((3,40))
ax1.set_ylabel(r'$\mathrm{\gamma\ cm^{-2}\ s^{-1}\ keV^{-1}}$')
ax1.legend()

for i in range(3):
	j = i+3
	ax2.errorbar(data[j][0],data[j][2], xerr=data[j][1], yerr=data[j][3], ls = 'none', color='C' + str(i))	#pl ratio
	ax2.axhline(y=1, lw = 0.8)
	# ax2.set_ylim((-0.8,0.8))
	ax2.set_yscale('log')
ax2.set_ylabel('data/model')
ax2.set_xlabel('Energy (keV)')
ax2.set_yticks([1, 100])

plt.tight_layout()
plt.savefig('/Users/sean/Desktop/HEAD2019/tbabs_diskbb_frozennH_xrt_nustar_30keV_ldata.pdf')
plt.close()


data = []
i = 0
with open('/Users/sean/Desktop/HEAD2019/nu90201034002A01_cl_tc_bc_regfilt_nustar_fpma_E3-40_pds.qdp', 'r') as file:
	for line in file:
		if i >0:
			temp = line.split()
			data.append(temp)
		i += 1

data = np.array(data).astype(float).T

fig, ax = plt.subplots(1, 1, figsize=(12,8))
ax.errorbar(data[0],data[1], yerr=data[2], ls = 'steps-mid')
# ax.ticklabel_format(axis='y', style='plain')
ax.set_xlim(0.003, 100)
ax.set_xscale('log')
# ax.set_yscale('log')
ax.set_ylabel('Leahy Power')
ax.set_xlabel('Frequency (Hz)')

plt.tight_layout()
plt.savefig('/Users/sean/Desktop/HEAD2019/nu90201034002A01_cl_tc_bc_regfilt_nustar_fpma_E3-40_pds.pdf')
plt.close()
