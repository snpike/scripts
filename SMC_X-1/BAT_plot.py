from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.gridspec import GridSpec

sns.set_context('paper', font_scale=1.2)
sns.set_style("whitegrid")
sns.set_palette("colorblind")

obs003 = [5.614510150675926E+04, 5.614545923824074E+04]
obs001 = [5.611328497898149E+04, 5.611431833546297E+04]

hdulist = fits.open('/Users/sean/Desktop/SMC_X-1_BAT_lc_2012_clean.fits')
plt.figure()
plt.errorbar(hdulist[1].data['TIME'], hdulist[1].data['RATE'], yerr = hdulist[1].data['ERROR'], fmt = 'none', lw = 0.5)
plt.ylabel('Counts ' + r'$\mathrm{cm}^{-2}$' + ' ' + r'$\mathrm{s}^{-1}$' + ' ' + '(15-50 keV)')
plt.xlabel('MJD')
plt.ylim(-0.005, 0.02)
plt.xlim(56000, 56250)
plt.axvspan(obs001[0], obs001[1], color='red')
plt.axvspan(obs003[0], obs003[1], color='red')
plt.savefig('/Users/sean/Desktop/SMC_X-1/letter/figures/SMC_X-1_BAT_lc_2012_clean.eps')