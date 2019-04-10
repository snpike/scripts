import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.gridspec import GridSpec
from scipy.optimize import curve_fit

sns.set_context('talk')
sns.set_style("ticks")
sns.set_palette("colorblind")

def superorbital(x, P, phi):
	return 0.09*np.square(np.sin((((x-57620)*2*np.pi)/P) + phi)) + 0.01

#distance (4Mpc) in cm
# distance = 1.234E25

#distance (2.9Mpc) in cm from Koribalski+ 2004, H=70km/s/Mpc, 
distance = 8.95E24
area = 4*np.pi*np.square(distance)

TSTART = np.array([56323.39520559, 57623.39701707, 57625.59070919, 57628.72335829, 57774.62891035, 57790.70553194, 57807.83040635, 58136.73824799, 58143.05138877, 58171.33118118, 58176.04965113, 58501.15146318, 58511.03708183, 58514.09455351, 58523.52816166, 58531.81129906])
TSTOP = np.array([56323.54015184, 57623.47477121, 57625.60727198, 57628.73035203, 57774.77127668, 57790.96502874, 57807.90840747, 58136.81154771, 58143.18028880, 58171.53446455, 58176.25111641, 58501.62219373, 58511.58062789, 58514.10970178, 58523.99649929, 58531.89301621])
TCENTER = (TSTART + TSTOP)/2

NUSTAR_START = np.array([56325.04493056, 56326.12565972, 56328.21215278, 57623.28882242])
NUSTAR_STOP = np.array([56325.42579890, 56326.97497714, 56328.97915509, 57624.33953571])
NUSTAR_CENTER = (NUSTAR_START + NUSTAR_STOP)/2

# # Counts/cm2/s for 2-10 keV
# flux = np.array([8.3452, 1.992, 1.9303, 1.174, 1.092, 0.7660, 1.774, 8.677, 5.196, 7.4698, 5.1014])
# uperr = np.array([0.3838, 0.191, 0.2337, 0, 0, 0, 0, 0, 0.401, 0.4912, 0.3886])
# lowerr = np.array([1.1142, 1.0446, 1.1793, 0, 0, 0, 0, 0, 1.353, 1.3928, 1.2514])

# nustar_flux = 3.7787
# nustar_uperr = 0.2033
# nustar_lowerr = 0.2907

# plt.figure()
# plt.errorbar(TCENTER, flux, xerr = (TSTOP-TSTART)/2, yerr = [uperr, lowerr], uplims= (uperr==0),fmt = '.')
# plt.errorbar(NUSTAR_CENTER, nustar_flux, xerr = (NUSTAR_STOP-NUSTAR_START)/2, yerr = [[nustar_uperr], [nustar_lowerr]], color='red', fmt='.')

# #plt.axvspan(NUSTAR_START, NUSTAR_STOP, color='red')
# plt.ylabel(r'$\mathrm{10^{-4}\,photon\,cm^{-2}\,s^{-1}\ (2-10\,keV)}$')
# plt.xlabel('MJD')
# plt.tight_layout()
# plt.savefig('/Users/sean/Desktop/CircinusULX5_lc_photonflux.pdf')

# erg/cm2/s for 2-10 keV
flux = np.array([5.0948, 1.042, 0.4869, 0.7916, 2.7334, 4.3208, 2.7502, 0.5487, 0.3771, 0.2774, 0.3110, 0.5885])
uperr = np.array([0.0292, 0, 0, 0, 0.0096, 0.0332, 0.00001, 0, 0, 0, 0, 0])
lowerr = np.array([2.4578, 0, 0, 0, 2.1814, 2.8828, 2.368, 0, 0, 0, 0, 0])

new_flux = np.sqrt((np.square(flux + uperr) + np.square(flux-lowerr))/2.0)
new_uperr = (flux + uperr)-new_flux
new_lowerr = new_flux - (flux-lowerr)

nustar_flux = np.array([5.271, 5.4978, 4.2581, 1.3558])
nustar_uperr = np.array([0.118, 0.0912, 0.0769, 0.0402])
nustar_lowerr = np.array([0.234, 0.1188, 0.1281, 0.0878])

ctrate = [1.051e-01, 2.774e-02, 1.478e-02, 4.297e-02, 5.628e-02, 9.084e-02, 6.683e-02, 1.932e-02, 1.124e-02, 8.614e-03, 1.362e-02, 2.375e-02, 1.666e-02, 1.785e-02, 2.768e-02, 8.298e-02]
ctrate_err = [5.495e-03, 3.865e-03, 4.977e-03, 9.161e-03, 4.648e-03, 6.465e-03, 5.330e-03, 3.278e-03, 2.504e-03, 2.349e-03, 2.894e-03, 2.928e-03, 3.180e-03, 3.848e-03, 3.147e-03, 5.773e-03]

lw = 1.7

popt, pcov = curve_fit(superorbital, TCENTER[1:], ctrate[1:], p0 = [75, 1.0])
print(popt)
print(pcov)

# plt.figure()
# plt.errorbar(TCENTER, new_flux, yerr = [new_lowerr, new_uperr], uplims= (uperr==0),fmt = '.', label = 'Swift XRT', linewidth = lw)
# plt.errorbar(NUSTAR_CENTER, nustar_flux, yerr = [nustar_lowerr, nustar_uperr], color='red', fmt='.', label = 'NuSTAR', linewidth = lw)

# plt.ylabel(r'$\mathrm{10^{-12}\,erg\,cm^{-2}\,s^{-1}\ (2-10\,keV)}$')
# plt.xlabel('MJD')
# plt.yscale('log')
# plt.legend()
# plt.tight_layout()
# plt.savefig('/Users/sean/Desktop/AstroData/Circinus_ULX5/figures/CircinusULX5_lc_ergflux.pdf')
# plt.close()

# plt.figure()
# plt.errorbar(TCENTER, new_flux*area*(10**(-12)), yerr = [new_lowerr*area*(10**(-12)), new_uperr*area*(10**(-12))], uplims= (uperr==0),fmt = '.', label = 'Swift XRT', linewidth=lw)
# plt.errorbar(NUSTAR_CENTER, nustar_flux*area*(10**(-12)), yerr = [nustar_lowerr*area*(10**(-12)), nustar_uperr*area*(10**(-12))], color='red', fmt='.', label='NuSTAR', linewidth=lw)

# plt.ylabel(r'$\mathrm{erg\,s^{-1}\ (2-10\,keV)}$')
# plt.yscale('log')
# plt.xlabel('MJD')
# plt.legend()
# plt.tight_layout()
# plt.savefig('/Users/sean/Desktop/AstroData/Circinus_ULX5/figures/CircinusULX5_lc_lum.pdf')
# plt.close()

plt.figure()
plt.errorbar(TCENTER[1:], ctrate[1:], yerr = ctrate_err[1:], fmt = '.', linewidth=lw)
plt.plot(np.arange(57620, 58600, 1), superorbital(np.arange(57620, 58600, 1), *popt))

plt.ylabel(r'XRT Count rate $(\mathrm{s^{-1}) (0.2-10\,keV)}$')
# plt.yscale('log')
plt.xlabel('MJD')
plt.tight_layout()
plt.show()
# plt.savefig('/Volumes/LaCie/AstroData/Circinus_ULX5/figures/CircinusULX5_XRT_ctrate.pdf')
plt.close()

# Tin_T = (np.array([56323.54015184, 57623.47477121, 57774.77127668, 57790.96502874, 57807.90840747, 58501.62219373]) + \
# 		np.array([56323.39520559, 57623.39701707, 57774.62891035, 57790.70553194, 57807.83040635, 58501.15146318]))/2.0

# Tin = [1.87, 1.34, 1.39, 1.69, 1.50, 0.72]

# Tin_up = [0.73, 1.20, 0.75, 0.74, 0.76, 0.39]

# Tin_low = [0.45, 0.50, 0.40, 0.41, 0.40, 0.22]

# Tin_T_nu = (np.array([56325.04493056, 56326.12565972, 56328.21215278, 57623.28882242]) + \
# 		np.array([56325.42579890, 56326.97497714, 56328.97915509, 57624.33953571]))/2.0

# Tin_nu = [1.8, 1.97, 1.72, 1.28]

# Tin_up_nu = [0.09, 0.06, 0.06, 0.07]

# Tin_low_nu = [0.08, 0.06, 0.06, 0.06]

# plt.figure()
# plt.errorbar(Tin_T, Tin, yerr = [Tin_low, Tin_up], fmt = '.', label = 'Swift XRT', linewidth=lw)
# plt.errorbar(Tin_T_nu, Tin_nu, yerr = [Tin_low_nu, Tin_up_nu], color='red', fmt='.', label='NuSTAR', linewidth=lw)

# plt.ylabel(r'$T_{in}\mathrm{\,(keV)}$')
# #plt.yscale('log')
# plt.xlabel('MJD')
# plt.legend(loc=3)
# plt.tight_layout()
# plt.savefig('/Users/sean/Desktop/AstroData/Circinus_ULX5/figures/CircinusULX5_Tin.pdf')
# plt.close()

# Tin = [0.64, 0.50, 0.49, 0.59, 0.58, 0.92]

# Tin_up = [0.26, 0.62, 0.35, 0.31, 0.38, 0.90]

# Tin_low = [0.19, 0.31, 0.22, 0.21, 0.23, 0.48]

# plt.figure()
# plt.errorbar(Tin_T, Tin, yerr = [Tin_low, Tin_up], fmt = '.', label = 'Swift XRT', linewidth=lw)

# plt.ylabel(r'$n_{\mathrm{H}}\,(10^{22}\,\mathrm{cm^{-2}})$')
# #plt.yscale('log')
# plt.xlabel('MJD')
# plt.legend(loc=3)
# plt.tight_layout()
# plt.savefig('/Users/sean/Desktop/AstroData/Circinus_ULX5/figures/CircinusULX5_nH.pdf')
# plt.close()
