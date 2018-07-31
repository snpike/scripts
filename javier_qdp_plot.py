import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.interpolate import interp1d
from sklearn.decomposition import PCA, NMF

sns.set_context('paper', font_scale=1.0, rc={"lines.linewidth": 0.8})
sns.set_style("ticks")
sns.set_palette("colorblind")

epochs = {'epochI': '/Users/sean/Desktop/for_sean/obs-1/obs-1_spec.qdp', \
		  'epochII': '/Users/sean/Desktop/for_sean/obs-2/obs-2_spec.qdp'}


file = open(epochs['epochI'], 'r')

data = [[]]
i = 0
j = 0
for line in file:
	if i >2:
		temp = line.split()
		if temp[0] != 'NO':
			data[j].append(temp)
		else:
			j += 1
			data.append([])
	i += 1

file.close()

data[0] = np.array(data[0]).astype(float).T
data[1] = np.array(data[1]).astype(float).T

spec1 = interp1d(data[0][0], data[0][2], kind = 'quadratic')

# plt.figure()
# plt.errorbar(data[0][0],data[0][2], xerr=data[0][1], yerr=data[0][3], ls = 'none', color='black', label='FPMA')	#spectrum
# plt.ylim((np.min(data[0][2]) * 0.8, np.max(data[0][2]) * 1.2))
# plt.errorbar(data[1][0],data[1][2], xerr=data[1][1], yerr=data[1][3], ls = 'none', color='red', label='FPMB')	#spectrum
# plt.plot(np.linspace(5, 50, num = 1000), spec1(np.linspace(5, 50, num = 1000)))
# plt.xscale('log')
# plt.yscale('log')
# plt.ylabel(r'$\mathrm{normalized counts\  s^{-1}\ keV^{-1}}$')
# plt.legend()

# plt.tight_layout(h_pad = 0)
# plt.show()
# plt.close()

file = open(epochs['epochII'], 'r')

data = [[]]
i = 0
j = 0
for line in file:
	if i >2:
		temp = line.split()
		if temp[0] != 'NO':
			data[j].append(temp)
		else:
			j += 1
			data.append([])
	i += 1

file.close()

data[0] = np.array(data[0]).astype(float).T
data[1] = np.array(data[1]).astype(float).T

spec2 = interp1d(data[0][0], data[0][2], kind = 'quadratic')

# plt.figure()
# plt.errorbar(data[0][0],data[0][2], xerr=data[0][1], yerr=data[0][3], ls = 'none', color='black', label='FPMA')	#spectrum
# plt.ylim((np.min(data[0][2]) * 0.8, np.max(data[0][2]) * 1.2))
# plt.errorbar(data[1][0],data[1][2], xerr=data[1][1], yerr=data[1][3], ls = 'none', color='red', label='FPMB')	#spectrum
# plt.plot(np.linspace(5, 50, num = 1000), spec2(np.linspace(5, 50, num = 1000)))
# plt.xscale('log')
# plt.yscale('log')
# plt.ylabel(r'$\mathrm{normalized counts\  s^{-1}\ keV^{-1}}$')
# plt.legend()

# plt.tight_layout(h_pad = 0)
# plt.show()
# plt.close()

X = np.array([spec1(np.linspace(5, 50, num = 1000)), spec2(np.linspace(5, 50, num = 1000))]).T
pca = PCA()
X_new = pca.fit_transform(X)
pca_comp = pca.components_
pca_sv = pca.singular_values_
pca_mean = pca.mean_

print('PCA evolution:')
print(pca_comp)

plt.figure()
for comp in (np.dot(X_new, pca_comp) + pca_mean).T:
	plt.plot(np.linspace(5, 50, num = 1000), comp, linewidth = 3.0)
for comp in X.T:
	plt.plot(np.linspace(5, 50, num = 1000), comp)
plt.xscale('log')
plt.yscale('log')
plt.show()
plt.close()

plt.figure()
for comp in (X_new + pca_mean).T:
	plt.plot(np.linspace(5, 50, num = 1000), comp, linewidth = 3.0)
plt.xscale('log')
plt.yscale('log')
plt.title('PCA components')
plt.show()
plt.close()

print(np.sum(np.dot(X_new, pca_comp) + pca_mean - X))

for i in range(10):
	nmf = NMF()
	W = nmf.fit_transform(X)
	H = nmf.components_
	# print(W)
	print('NMF evolution:')
	print(H)


	plt.figure()
	for comp in W.T:
		plt.plot(np.linspace(5, 50, num = 1000), comp)
	plt.xscale('log')
	plt.yscale('log')
	plt.title('NMF components')
	plt.show()
	plt.close()

	# plt.figure()
	# for comp in H:
	# 	plt.plot(range(len(comp)), comp)
	# plt.show()
	# plt.close()

	plt.figure()
	for comp in np.dot(W, H).T:
		plt.plot(np.linspace(5, 50, num = 1000), comp)
	# plt.plot(np.linspace(5, 50, num = 1000), spec1(np.linspace(5, 50, num = 1000)))
	# plt.plot(np.linspace(5, 50, num = 1000), spec2(np.linspace(5, 50, num = 1000)))
	plt.xscale('log')
	plt.yscale('log')
	plt.show()
	plt.close()

	print()

