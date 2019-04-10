from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.gridspec import GridSpec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, InsetPosition, mark_inset

sns.set_context('paper', font_scale=1.5)
#sns.set_context('poster')
sns.set_style("ticks")
sns.set_palette("colorblind")

def moving_average(x, n=5):
	# Determine the running average of y(x) using stepsize (in indices, not units of x)
    new_x = np.cumsum(x, dtype=float)
    new_x[n:] = new_x[n:] - new_x[:-n]
    return new_x[n-1:]/n


obs003 = [5.614510150675926E+04, 5.614545923824074E+04]
obs001 = [5.611328497898149E+04, 5.611431833546297E+04]

hdulist = fits.open('/Users/sean/AstroData/SMC_X-1_BAT_lc_2012_clean.fits')
fig, ax1 = plt.subplots(figsize=(8, 4.5))

# print(np.argwhere(hdulist[1].data['TIME']>56080))
# print(np.argwhere(hdulist[1].data['TIME']<56180))
ax1.errorbar(hdulist[1].data['TIME']-56000, hdulist[1].data['RATE'], yerr = hdulist[1].data['ERROR'], fmt = 'none', lw = 0.5)
#ax1.subplots_adjust(left=0.15)
ax1.set_ylabel('Counts ' + r'$\mathrm{cm}^{-2}$' + ' ' + r'$\mathrm{s}^{-1}$' + ' ' + '(15-50 keV)')
ax1.set_xlabel('MJD - 56000')
ax1.set_ylim(-0.005, 0.02)
ax1.set_xlim(0, 250)
ax1.axvspan(obs001[0]-56000, obs001[1]-56000, color='red')
ax1.axvspan(obs003[0]-56000, obs003[1]-56000, color='red')
ax1.plot(moving_average(hdulist[1].data['TIME'], n=8)-56000, moving_average(hdulist[1].data['RATE'], n=8), lw=2.0)
ax1.tick_params(direction='in')
plt.tight_layout()

# ax2 = plt.axes([0,0,1,1])
# # Manually set the position and relative size of the inset axes within ax1
# ip = InsetPosition(ax1, [0.08,0.1,0.3,0.3])
# ax2.set_axes_locator(ip)
# ax2.plot(moving_average(hdulist[1].data['TIME'][150:250], n=8)-56000, moving_average(hdulist[1].data['RATE'][150:250], n=8))
# ax2.axvspan(obs001[0]-56000, obs001[1]-56000, color='red')
# ax2.axvspan(obs003[0]-56000, obs003[1]-56000, color='red')

#plt.show()
#plt.savefig('/Users/sean/Desktop/SMC_X-1/letter/figures/SMC_X-1_BAT_lc_2012_clean.eps')
# plt.tight_layout()
# plt.show()
plt.savefig('/Users/sean/Desktop/SMC_X-1_BAT_lc_2012_clean.eps')
plt.close()