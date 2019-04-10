import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from tqdm import tqdm
import scipy.stats as stats

def line(x, m, b):
	return (m*x) + b

def line_inv(y, m, b):
	return (y-b)/m


reg_width = 109

src_cntrt = 6.189E-3
#src_cntrt = 1.238E-1

src_sigma = 7

bkg_cntrt = 0.02

exposure = [50*(i+1) for i in range(6)] + [400 + (100*i) for i in range(12)]
print(exposure)

false = [0 for i in range(18)]

ntrials = 10000
one_sigma = 0.317310507863
three_sigma = 0.002699796063
four_sigma = 0.00006334
five_sigma = 0.000000573303

for i in range(len(exposure)):
	print(i+1)
	t = exposure[i]
	src_counts_array = np.random.poisson(src_cntrt*t, ntrials)
	#expected_A_array = src_counts_array/(src_sigma*np.sqrt(2*np.pi))
	for j in tqdm(range(ntrials)):
		src_counts = src_counts_array[j]
		#expected_A = expected_A_array[j]

		region = np.random.poisson(bkg_cntrt*t/(reg_width*reg_width), (reg_width, reg_width))
		y, x = np.mgrid[:reg_width, :reg_width]
		# plt.figure()
		# plt.imshow(region)
		# plt.show()

		p_init = models.Gaussian2D()
		fit_p = fitting.LevMarLSQFitter()
		p = fit_p(p_init, x, y, region)

		if (p.x_stddev < reg_width) and (p.y_stddev < reg_width):
			expected_A = 2*src_counts/((p.x_stddev + p.y_stddev)*(np.sqrt(2*np.pi)))
			if (p.amplitude >= expected_A):
				false[i] += 1
false_rate = np.array(false)/ntrials
print(false_rate)
exposure = np.array(exposure)

#false rate for low state
# false_rate = [0.2804, 0.0807, 0.0293, 0.0153, 0.0108, 0.0087, 0.0083, 0.0084, 0.0071, 0.0086,
# 0.006,  0.0062, 0.0052, 0.0049, 0.0044, 0.0045, 0.0053, 0.0036]

# false rate for high state
# false_rate = [1.08e-03, 8.00e-05, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 2.00e-05, 2.00e-05,
#  2.00e-05, 2.00e-05, 5.00e-05, 6.00e-05, 7.00e-05, 1.20e-04, 9.00e-05, 4.00e-05,
#  1.20e-04, 1.50e-04]

coeff = np.polyfit(exposure, false_rate, 20)
polyf = np.poly1d(coeff)
print(polyf.c)
print(polyf.o)

plt.figure()
plt.plot(exposure, false_rate)
#plt.plot(exposure, polyf(exposure))
plt.axhline(three_sigma, color='red')
plt.axhline(four_sigma, color='yellow')
plt.axhline(five_sigma, color='green')
plt.show()
plt.close()

# onesig_coeff = coeff
# onesig_coeff[-1] -= one_sigma
# onesig_poly = np.poly1d(onesig_coeff)
# print('One sigma detection at: ' + str(np.roots(onesig_coeff)))

threesig_coeff = coeff
threesig_coeff[-1] += one_sigma
threesig_coeff[-1] -= three_sigma
threesig_poly = np.poly1d(threesig_coeff)
print('Three sigma detection at: ' + str(np.roots(threesig_coeff)))

