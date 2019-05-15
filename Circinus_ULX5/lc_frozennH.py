import csv
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import odr
import seaborn as sns
import scipy.signal as sig

sns.set_context('poster')
sns.set_style('ticks')
sns.set_palette('colorblind')

#distance (3Mpc) in cm from Koribalski+ 2004, H=70km/s/Mpc, 
distance = 9.3E24
area = 4*np.pi*np.square(distance)

def superorbital(x, P, phi):
	return 0.06*np.square(np.cos((np.pi*x/P) + phi)) + 0.01

# def superorbital(x, P, phi, C):
# 	return A*np.square(np.cos((np.pi*x/P) + phi)) + C

def powerlaw(beta, x):
	A, n = beta
	return A*np.power(x, n)

def powerlaw_inv(beta, y):
	A, n = beta
	return np.power(y/A, 1/n)


data = {}

with open('/Users/sean/Desktop/CircinusULX5_xrt/diskbb, nH = 0.6 frozen (1-sigma err)-s.csv', 'r') as file:
	reader = csv.DictReader(file)
	for row in reader:
		for col in row:
			if col not in data:
				data[col] = []
			data[col].append(row[col])

# print(data)

flux = np.array([float(x) for x in data['Flux (erg/cm^2/s) (2-10 keV)']])
flux_p = np.array([float(x) for x in data['Flux uncertainty +']])
flux_m = np.array([float(x) for x in data['Flux uncertainty -']])
TSTART = np.array([float(x) for x in data['TSTART (MJD)']])
TSTOP = np.array([float(x) for x in data['TSTOP (MJD)']])
T = np.array([float(x) for x in data['Tin (keV)']])
T_p = np.array([float(x) for x in data['Tin uncertainty +']])
T_m = np.array([float(x) for x in data['Tin uncertainty -']])
TCENTER = (TSTOP + TSTART)/2

NUSTAR_mask = np.nonzero(np.array(data['Instrument']) == 'NuSTAR')
XRT_mask = np.nonzero(np.array(data['Instrument']) == 'Swift')
fluxlim_mask = np.nonzero(flux[XRT_mask] == flux_m[XRT_mask])
fluxmeasure_mask = np.nonzero(flux[XRT_mask] != flux_m[XRT_mask])
T_mask = np.nonzero(T_p != 999)
T_nustar = np.nonzero((T_p != 999) * (np.array(data['Instrument']) == 'NuSTAR'))
T_xrt = np.nonzero((T_p != 999) * (np.array(data['Instrument']) == 'Swift'))
Tandflux = np.nonzero((T_p != 999) * (flux != flux_m))
Tandflux_xrt = np.nonzero((T_p != 999) * (flux != flux_m) * (np.array(data['Instrument']) == 'Swift'))
Tandflux_nustar = np.nonzero((T_p != 999) * (flux != flux_m) * (np.array(data['Instrument']) == 'NuSTAR'))

# new_flux_p = flux_p
# new_flux_p[fluxlim_mask] = 0.0
# new_flux_m = flux_m
# new_flux_m[fluxlim_mask] = 0.0

fig, ax1 = plt.subplots(1,1, figsize=(12, 8))
ax1.errorbar(TCENTER[NUSTAR_mask], flux[NUSTAR_mask], yerr= [flux_m[NUSTAR_mask], flux_p[NUSTAR_mask]], fmt = '.', color = 'C1', label = 'NuSTAR')
ax1.errorbar(TCENTER[XRT_mask][fluxmeasure_mask], flux[XRT_mask][fluxmeasure_mask], yerr= [flux_m[XRT_mask][fluxmeasure_mask], flux_p[XRT_mask][fluxmeasure_mask]], fmt = '.', color = 'C0', label = 'Swift XRT')
ax1.errorbar(TCENTER[XRT_mask][fluxlim_mask], flux[XRT_mask][fluxlim_mask], yerr = 0, uplims = flux[XRT_mask][fluxlim_mask], fmt = '.', color = 'C0')
ax1.legend()
# plt.xlim((57550, 58600))
# plt.ylim(1.3e-13, 7e-12)
ax1.set_ylabel(r'2-10 keV Flux ($erg\,cm^2\,s^{-1}$)')
ax1.set_xlabel('MJD')
ax1.set_yscale('log')
ax2 = ax1.twinx()
ax2.set_ylim(ax1.get_ylim()[0]*area, ax1.get_ylim()[1]*area)
ax2.set_ylabel(r'2-10 keV Luminosity ($erg\,s^{-1}$)')
ax2.set_yscale('log')
plt.tight_layout()
# plt.savefig('/Users/sean/Desktop/HEAD2019/flux.pdf')
plt.show()
plt.close()


fig, ax1 = plt.subplots(1,1, figsize=(12, 8))
ax1.errorbar(TCENTER[NUSTAR_mask], flux[NUSTAR_mask], yerr= [flux_m[NUSTAR_mask], flux_p[NUSTAR_mask]], fmt = '.', color = 'C1', label = 'NuSTAR')
ax1.errorbar(TCENTER[XRT_mask][fluxmeasure_mask], flux[XRT_mask][fluxmeasure_mask], yerr= [flux_m[XRT_mask][fluxmeasure_mask], flux_p[XRT_mask][fluxmeasure_mask]], fmt = '.', color = 'C0', label = 'Swift XRT')
ax1.errorbar(TCENTER[XRT_mask][fluxlim_mask], flux[XRT_mask][fluxlim_mask], yerr = 0, uplims = flux[XRT_mask][fluxlim_mask], fmt = '.', color = 'C0')
ax1.legend()
ax1.set_yscale('log')
ax1.set_xlim((57550, 58600))
ax1.set_ylim(1.3e-13, 7e-12)
ax1.set_ylabel(r'2-10 keV Flux ($erg\,cm^2\,s^{-1}$)')
ax1.set_xlabel('MJD')
ax2 = ax1.twinx()
ax2.set_ylim(ax1.get_ylim()[0]*area, ax1.get_ylim()[1]*area)
ax2.set_ylabel(r'2-10 keV Luminosity ($erg\,s^{-1}$)')
ax2.set_yscale('log')
plt.tight_layout()
# plt.savefig('/Users/sean/Desktop/HEAD2019/flux_new.pdf')
plt.show()
plt.close()

plt.figure(figsize=(12,8))
plt.errorbar(TCENTER[NUSTAR_mask], flux[NUSTAR_mask], yerr= [flux_m[NUSTAR_mask], flux_p[NUSTAR_mask]], fmt = '.', color = 'C1', label = 'NuSTAR')
plt.errorbar(TCENTER[XRT_mask][fluxmeasure_mask], flux[XRT_mask][fluxmeasure_mask], yerr= [flux_m[XRT_mask][fluxmeasure_mask], flux_p[XRT_mask][fluxmeasure_mask]], fmt = '.', color = 'C0', label = 'Swift XRT')
plt.errorbar(TCENTER[XRT_mask][fluxlim_mask], flux[XRT_mask][fluxlim_mask], yerr = 0, uplims = flux[XRT_mask][fluxlim_mask], fmt = '.', color = 'C0')
plt.legend()
plt.ylabel(r'2-10 keV Flux ($erg\,cm^2\,s^{-1}$)')
plt.xlabel('MJD')
plt.xlim((57600, 57830))
plt.ylim(1.2e-13, 6e-12)
plt.yscale('log')
plt.tight_layout()
plt.show()
plt.close()

plt.figure(figsize=(12,8))
plt.errorbar(TCENTER[NUSTAR_mask], flux[NUSTAR_mask], yerr= [flux_m[NUSTAR_mask], flux_p[NUSTAR_mask]], fmt = '.', color = 'C1', label = 'NuSTAR')
plt.errorbar(TCENTER[XRT_mask][fluxmeasure_mask], flux[XRT_mask][fluxmeasure_mask], yerr= [flux_m[XRT_mask][fluxmeasure_mask], flux_p[XRT_mask][fluxmeasure_mask]], fmt = '.', color = 'C0', label = 'Swift XRT')
plt.errorbar(TCENTER[XRT_mask][fluxlim_mask], flux[XRT_mask][fluxlim_mask], yerr = 0, uplims = flux[XRT_mask][fluxlim_mask], fmt = '.', color = 'C0')
plt.legend()
plt.ylabel(r'2-10 keV Flux ($erg\,cm^2\,s^{-1}$)')
plt.xlabel('MJD')
plt.xlim((58130, 58180))
plt.ylim(1.2e-13, 6e-12)
plt.yscale('log')
plt.tight_layout()
plt.show()
plt.close()

plt.figure(figsize=(12,8))
plt.errorbar(TCENTER[NUSTAR_mask], flux[NUSTAR_mask], yerr= [flux_m[NUSTAR_mask], flux_p[NUSTAR_mask]], fmt = '.', color = 'C1', label = 'NuSTAR')
plt.errorbar(TCENTER[XRT_mask][fluxmeasure_mask], flux[XRT_mask][fluxmeasure_mask], yerr= [flux_m[XRT_mask][fluxmeasure_mask], flux_p[XRT_mask][fluxmeasure_mask]], fmt = '.', color = 'C0', label = 'Swift XRT')
plt.errorbar(TCENTER[XRT_mask][fluxlim_mask], flux[XRT_mask][fluxlim_mask], yerr = 0, uplims = flux[XRT_mask][fluxlim_mask], fmt = '.', color = 'C0')
# plt.legend()
plt.ylabel(r'2-10 keV Flux ($erg\,cm^2\,s^{-1}$)')
plt.xlabel('MJD')
plt.xlim((58495, 58560))
plt.ylim(1.2e-13, 6e-12)
plt.yscale('log')
plt.tight_layout()
plt.show()
# plt.savefig('/Users/sean/Desktop/HEAD2019/flux_2019.pdf')
plt.close()

xrt_cntrt = np.array([float(x) for x in data['XRT Count Rate (0.3-10 keV)'][5:]])
xrt_cntrt_uncrt = np.array([float(x) for x in data['XRT Count Rate uncertainty +/-'][5:]])
plt.figure(figsize=(12,8))
plt.errorbar(TCENTER[5:], xrt_cntrt, yerr= xrt_cntrt_uncrt, fmt = '.')
plt.yscale('log')
plt.xlabel('MJD')
plt.ylabel(r'Swift XRT 0.3-10 keV count rate ($ct\,s^{-1}$)')
plt.tight_layout()
plt.show()
# plt.savefig('/Users/sean/Desktop/HEAD2019/xrtrate_new.pdf')
plt.close()

plt.figure(figsize=(12,8))
plt.errorbar(TCENTER[5:], xrt_cntrt, yerr= xrt_cntrt_uncrt, fmt = '.', color = 'C1')
plt.xlim((58495, 58560))
plt.tight_layout()
plt.show()
plt.close()


realdata = odr.RealData(T[Tandflux], flux[Tandflux], sx=np.sqrt((np.square(T_m[Tandflux]) + np.square(T_p[Tandflux]))/2), sy=np.sqrt((np.square(flux_m[Tandflux]) + np.square(flux_p[Tandflux]))/2))
model = odr.Model(powerlaw)
myodr = odr.ODR(realdata, model, beta0=[1e-12, 4])
myoutput = myodr.run()
myoutput.pprint()

plt.figure(figsize=(12,8))
plt.errorbar(T[Tandflux_nustar], flux[Tandflux_nustar], xerr=[T_m[Tandflux_nustar], T_p[Tandflux_nustar]], yerr=[flux_m[Tandflux_nustar], flux_p[Tandflux_nustar]], fmt='.', label='NuSTAR', color='C1')
plt.errorbar(T[Tandflux_xrt], flux[Tandflux_xrt], xerr=[T_m[Tandflux_xrt], T_p[Tandflux_xrt]], yerr=[flux_m[Tandflux_xrt], flux_p[Tandflux_xrt]], fmt='.', label='Swift XRT', color='C0')
plt.plot(np.arange(0.6, 2.5, 0.1), powerlaw(myoutput.beta, np.arange(0.6, 2.5, 0.1)), color='C2', label = r'$F \propto T^{2.93 \pm 0.24}$')
plt.legend()
plt.ylabel(r'2-10 keV Flux ($erg\,cm^2\,s^{-1}$)')
plt.xlabel(r'$kT_{in} (keV)$')
plt.yscale('log')
plt.xscale('log')
plt.tight_layout()
plt.show()
# plt.savefig('/Users/sean/Desktop/HEAD2019/Tin_flux.pdf')
plt.close()

print(powerlaw_inv(myoutput.beta, 1e-12))

plt.figure(figsize=(12,8))
plt.errorbar(TCENTER[T_nustar], T[T_nustar], yerr=[T_m[T_nustar], T_p[T_nustar]], fmt='.', label='NuSTAR', color='C1')
plt.errorbar(TCENTER[T_xrt], T[T_xrt], yerr=[T_m[T_xrt], T_p[T_xrt]], fmt='.', label='Swift XRT', color='C0')
plt.legend()
plt.xlabel('MJD')
plt.ylabel(r'$kT_{in} (keV)$')
plt.tight_layout()
plt.show()
# plt.savefig('/Users/sean/Desktop/HEAD2019/Tin_MJD.pdf')
plt.close()

plt.figure(figsize=(12,8))
plt.errorbar(TCENTER[T_nustar][3:], T[T_nustar][3:], yerr=[T_m[T_nustar][3:], T_p[T_nustar][3:]], fmt='.', label='NuSTAR', color='C1')
plt.errorbar(TCENTER[T_xrt][1:], T[T_xrt][1:], yerr=[T_m[T_xrt][1:], T_p[T_xrt][1:]], fmt='.', label='Swift XRT', color='C0')
plt.legend()
plt.xlabel('MJD')
plt.ylabel(r'$kT_{in} (keV)$')
plt.tight_layout()
plt.show()
# plt.savefig('/Users/sean/Desktop/HEAD2019/Tin_MJD_new.pdf')
plt.close()

# popt, pcov = curve_fit(superorbital, TCENTER[-6:]-np.mean(TCENTER[-6:]), xrt_cntrt[-6:], sigma = xrt_cntrt_uncrt[-6:], p0 = [0.06, 60, np.pi/3, 0.01], bounds=([0, 50, 0, 0], [1.0, 70, np.pi, 0.1]))
popt, pcov = curve_fit(superorbital, TCENTER[15:]-np.mean(TCENTER[15:]), xrt_cntrt[10:], sigma = xrt_cntrt_uncrt[10:], p0 = [60, np.pi/2], bounds=([10, 0], [100, np.pi]))
print(popt)
print(pcov)
print(np.diag(pcov))

# plt.figure(figsize=(12,8))
plt.figure(figsize=(16,8))
plt.errorbar(TCENTER[15:], xrt_cntrt[10:], yerr= xrt_cntrt_uncrt[10:], fmt = '.', label = 'Observed count rate')
plt.plot(np.arange(TCENTER[15], TCENTER[-1] + 0.5, 0.5), superorbital(np.arange(TCENTER[15]-np.mean(TCENTER[15:]), TCENTER[-1]-np.mean(TCENTER[15:]) + 0.5, 0.5), *popt), label = 'Model (P=' + str(int(round(popt[0], 0))) + ' days)', color='C2')
# plt.plot([58500, 58555], [0.05,0.05], linestyle='--', label = r'Trigger threshold (0.05 $\mathrm{ct\,s^{-1}}$)')
plt.legend()
plt.xlabel('MJD')
plt.ylabel(r'Swift XRT 0.3-10 keV count rate ($ct\,s^{-1}$)')
# plt.yscale('log')
plt.tight_layout()
# plt.show()
plt.savefig('/Users/sean/Desktop/CircinusULX5/xrtrate_2019.pdf')
plt.close()

plt.hist(xrt_cntrt, bins = 20)
plt.show()
plt.close()

# popt, pcov = curve_fit(superorbital, TCENTER[-16:-10], xrt_cntrt[-16:-10], sigma = xrt_cntrt_uncrt[-16:-10], p0 = [0.06, 65, np.pi/3, 0.01], bounds=([0, 10, 0, 0], [1.0, 100, 2*np.pi, 0.1]))
# print(popt)
# print(np.diag(pcov))

# plt.figure(figsize=(12,8))
# plt.errorbar(TCENTER[-16:-10], xrt_cntrt[-16:-10], yerr= xrt_cntrt_uncrt[-16:-10], fmt = '.')
# plt.plot(np.arange(TCENTER[-16], TCENTER[-11] + 0.5, 0.5), superorbital(np.arange(TCENTER[-16], TCENTER[-11] + 0.5, 0.5), *popt))
# plt.tight_layout()
# plt.show()
# plt.close()

# See how often we fit to 60 day period:
# ntrials = 10000
# popt_sim = []
# pcov_sim = []
# Pseeds = np.random.uniform(low=40,high=80, size=ntrials)
# # phiseeds = np.random.uniform(low=0,high=np.pi, size=ntrials)
# for n in range(ntrials):
# 	popt, pcov = curve_fit(superorbital, TCENTER[-6:]-np.mean(TCENTER[-6:]), xrt_cntrt[-6:], sigma = xrt_cntrt_uncrt[-6:], p0 = [Pseeds[n], np.pi/2], bounds=([10, 0], [100, np.pi]))
# 	popt_sim.append(popt)
# 	pcov_sim.append(np.diag(pcov))

# popt_sim = np.array(popt_sim)
# pcov_sim = np.array(pcov_sim)
# plt.figure(figsize=(12,8))
# plt.hist(popt_sim.T[0])
# plt.show()
# plt.close()

# plt.figure(figsize=(12,8))
# plt.scatter(popt_sim.T[0], pcov_sim.T[0])
# plt.show()
# plt.close()

# plt.figure(figsize=(12,8))
# plt.scatter(Pseeds, popt_sim.T[0])
# plt.show()
# plt.close()

# Simulate a light curve with superorbital modulation and see how many observations are needed:
# monitoring = [50,100,150,200]
# for l in monitoring:
# 	ntrials = 1000
# 	popt_sim = []
# 	pcov_sim = []
# 	rms = []
# 	print(l)

# 	for n in range(ntrials):
# 		sim_time=np.arange(0,l,10)
# 		sim_counts = (superorbital(sim_time, 61,  2.13)*3000)
# 		sim_lc = np.random.poisson(sim_counts)/3000
# 		popt, pcov = curve_fit(superorbital, sim_time, sim_lc, sigma = np.sqrt(sim_counts)/3000, p0 = [np.random.uniform(low=50, high=70), np.pi/2], bounds=([10, 0], [100, np.pi]))
# 		popt_sim.append(popt)
# 		pcov_sim.append(np.diag(pcov))
# 		rms.append(np.sqrt(np.sum(np.square(sim_lc))/len(sim_lc)))

# 	popt_sim = np.array(popt_sim)
# 	pcov_sim = np.array(pcov_sim)

# 	plt.figure(figsize=(12,8))
# 	plt.hist(rms)
# 	plt.show()
# 	plt.close()


# 	plt.figure(figsize=(12,8))
# 	plt.hist(popt_sim.T[0])
# 	plt.show()
# 	plt.close()

# 	plt.figure(figsize=(12,8))
# 	plt.hist(pcov_sim.T[0])
# 	plt.show()
# 	plt.close()

# 	# print(popt)
# 	# print(np.diag(pcov))
# 	# plt.figure(figsize=(12,8))
# 	# plt.errorbar(sim_time, sim_lc, yerr = np.sqrt((sim_lc)/3000), fmt='.')
# 	# plt.plot(np.arange(0, l + 0.5, 0.5), superorbital(np.arange(0, l + 0.5, 0.5), *popt))
# 	# plt.show()
# 	# plt.close()

# 	freqs = np.arange(0.001, 0.1, 0.001)
# 	lomb = sig.lombscargle(sim_time, sim_lc, freqs)
# 	plt.figure(figsize=(12,8))
# 	plt.plot(freqs, lomb)
# 	plt.show()
# 	plt.close()

# monitoring = [50,100,150,200]
# for l in monitoring:
# 	ntrials = 1000
# 	popt_sim = []
# 	pcov_sim = []
# 	rms= []
# 	print(l)

# 	for n in range(ntrials):
# 		sim_time=np.arange(0,l,10)
# 		sim_counts = np.random.uniform(low=40, high=200, size = len(sim_time))
# 		sim_lc = np.random.poisson(sim_counts)/3000
# 		rms.append(np.sqrt(np.sum(np.square(sim_lc))/len(sim_lc)))
# 	# 	popt, pcov = curve_fit(superorbital, sim_time, sim_lc, sigma = np.sqrt(sim_counts)/3000, p0 = [np.random.uniform(low=50, high=70), np.pi/2], bounds=([10, 0], [100, np.pi]))
# 	# 	popt_sim.append(popt)
# 	# 	pcov_sim.append(np.diag(pcov))

# 	plt.figure(figsize=(12,8))
# 	plt.hist(rms)
# 	plt.show()
# 	plt.close()

# 	# popt_sim = np.array(popt_sim)
# 	# pcov_sim = np.array(pcov_sim)

# 	# plt.figure(figsize=(12,8))
# 	# plt.hist(popt_sim.T[0])
# 	# plt.show()
# 	# plt.close()

# 	# plt.figure(figsize=(12,8))
# 	# plt.hist(pcov_sim.T[0])
# 	# plt.show()
# 	# plt.close()

# 	# print(popt)
# 	# print(np.diag(pcov))
# 	# plt.figure(figsize=(12,8))
# 	# plt.errorbar(sim_time, sim_lc, yerr = np.sqrt((sim_lc)/3000), fmt='.')
# 	# plt.plot(np.arange(0, l + 0.5, 0.5), superorbital(np.arange(0, l + 0.5, 0.5), *popt))
# 	# plt.show()
# 	# plt.close()

# 	freqs = np.arange(0.001, 0.1, 0.001)
# 	lomb = sig.lombscargle(sim_time, sim_lc, freqs)
# 	plt.figure(figsize=(12,8))
# 	plt.plot(freqs, lomb)
# 	plt.show()
# 	plt.close()


# sim_time=np.arange(0,10000,10)
# sim_counts = (superorbital(sim_time, 61,  2.13)*3000)
# sim_lc = np.random.poisson(sim_counts)/3000
# plt.figure(figsize=(12,8))
# plt.hist(sim_lc)
# plt.show()
# plt.close()


# sim_time=np.arange(0,10000,10)
# sim_counts = np.random.uniform(low=40, high=200, size = len(sim_time))
# sim_lc = np.random.poisson(sim_counts)/3000
# plt.figure(figsize=(12,8))
# plt.hist(sim_lc)
# plt.show()
# plt.close()

# popt, pcov = curve_fit(superorbital, TCENTER[-10:-6], xrt_cntrt[-10:-6], sigma = xrt_cntrt_uncrt[-10:-6], p0 = [0.06, 60, 0, 0.01], bounds=([0, 10, 0, 0], [1.0, 100, 2*np.pi, 0.1]))
# print(popt)
# print(np.diag(pcov))

# plt.figure()
# plt.errorbar(TCENTER[-10:-6], xrt_cntrt[-10:-6], yerr= xrt_cntrt_uncrt[-10:-6], fmt = '.', color = 'blue')
# plt.plot(np.arange(TCENTER[-10], TCENTER[-7] + 0.5, 0.5), superorbital(np.arange(TCENTER[-10], TCENTER[-7] + 0.5, 0.5), *popt))
# plt.show()
# plt.close()





