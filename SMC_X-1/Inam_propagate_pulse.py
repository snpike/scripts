import matplotlib.pyplot as plt
import numpy as np
import scipy.stats
import scipy.special
import scipy.fftpack
from scipy.optimize import curve_fit
from scipy.odr import ODR, RealData, Model

MJD_m50000 = []
freq = []
freq_err = []

def f(B, x):
	return B[0]*x + B[1]

with open('/Users/sean/Desktop/Inam10_SMCX1_pulsehistory.txt', 'r') as file:
	for line in file:
		splitline = line.split()
		# MJD is pretty straightforward, we just need to combine the digits
		MJD_m50000.append(float(splitline[1] + splitline[2])-50000)
		
		# frequency takes some more thought. Split the error from the measurement.
		# Then we need to determine how many sigfigs
		temp_freq = splitline[3][:-1].split('(')
		freq.append(float(temp_freq[0]))
		err_power = 2-len(temp_freq[0])
		freq_err.append(float(temp_freq[1])*(10**err_power))

# print(MJD_m50000)
# print('\n')
# print(freq)
# print('\n')
# print(freq_err)
# print('\n')

data = RealData(MJD_m50000, freq, sy = freq_err)
linear = Model(f)

# initial guess from Inam 2010 fit
myodr = ODR(data, linear, beta0 = [2.29257057E-6, 1.413])
myoutput = myodr.run()
myoutput.pprint()
print('\n')

Tref_III_m50000 = 56145.10372569-50000

print(myoutput.beta[0]/(24*3600))
print(myoutput.sd_beta[0]*1.644854/(24*3600))
print('\n')
print(1.0/f(myoutput.beta, Tref_III_m50000))
print(np.sqrt(np.square(Tref_III_m50000*myoutput.sd_beta[0]) + np.square(myoutput.sd_beta[1]))/np.square(f(myoutput.beta, Tref_III_m50000)))
