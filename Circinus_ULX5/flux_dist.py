import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sig
from scipy.optimize import curve_fit

def superorbital(x, P, phi):
	lc = 0.1 * np.square(np.cos((2*np.pi*x/P) + phi)) + 0.005
	return np.random.normal(lc, 0.01)

# def propellor(x, min, max):
# 	tmin = 30
# 	lc = []
# 	t = 0
# 	for xi in x:
# 		while t < np.normal(tmin, 5):


TSTART = np.array([56323.39520559, 57623.39701707, 57625.59070919, 57628.72335829, 57774.62891035, 57790.70553194, 57807.83040635, 58136.73824799, 58143.05138877, 58171.33118118, 58176.04965113, 58501.15146318, 58511.03708183, 58514.09455351, 58523.52816166, 58531.81129906])
TSTOP = np.array([56323.54015184, 57623.47477121, 57625.60727198, 57628.73035203, 57774.77127668, 57790.96502874, 57807.90840747, 58136.81154771, 58143.18028880, 58171.53446455, 58176.25111641, 58501.62219373, 58511.58062789, 58514.10970178, 58523.99649929, 58531.89301621])
TCENTER = (TSTART + TSTOP)/2

ctrate = [1.051e-01, 2.774e-02, 1.478e-02, 4.297e-02, 5.628e-02, 9.084e-02, 6.683e-02, 1.932e-02, 1.124e-02, 8.614e-03, 1.362e-02, 2.375e-02, 1.666e-02, 1.785e-02, 2.768e-02, 8.298e-02]
ctrate_err = [5.495e-03, 3.865e-03, 4.977e-03, 9.161e-03, 4.648e-03, 6.465e-03, 5.330e-03, 3.278e-03, 2.504e-03, 2.349e-03, 2.894e-03, 2.928e-03, 3.180e-03, 3.848e-03, 3.147e-03, 5.773e-03]


P = np.arange(20, 100, 1)
freq = np.arange(0.01, 0.05, 0.0001)

sim_sup = superorbital(np.arange(0,1000,0.1), 75, 0)
# sim_prop = 


lw = 1.7
plt.figure()
plt.errorbar(TCENTER[-5:], ctrate[-5:], yerr = ctrate_err[-5:], fmt = '.', linewidth=lw)
# plt.plot(np.arange(0,1000,0.1), sim_sup)
plt.ylabel(r'XRT Count rate $(\mathrm{s^{-1}) (0.2-10\,keV)}$')
# plt.yscale('log')
plt.xlabel('MJD')
plt.tight_layout()
plt.show()
plt.close()

plt.figure()
superorb_dist = plt.hist(sim_sup)
plt.show()
plt.close()

plt.figure()
plt.hist(ctrate)
plt.show()
plt.close()

pgram = sig.lombscargle(TCENTER, ctrate, freq, precenter=True)

plt.figure()
plt.plot(1/freq, pgram)
plt.show()
plt.close()