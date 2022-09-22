import matplotlib.pyplot as plt 


atoll = [1e-2, 5e-1]
Z = [5e-1, 3.0]
be = [1e-4, 1.0] 
sgxb = [1e-2, 5.0]
ulxp = [10.,500.]


plt.figure(figsize=(9,4))
plt.plot(atoll, [5.0,5.0], marker='|', lw = 5.0)
plt.plot(Z, [4.0,4.0], marker='|', lw = 5.0)
plt.plot(be, [3.0,3.0], marker='|', lw = 5.0)
plt.plot(sgxb, [2.0,2.0], marker='|', lw = 5.0)
plt.plot(ulxp, [1.0,1.0], marker='|', lw = 5.0)
plt.axhspan(3.5, 6, facecolor='blue', alpha=0.2,zorder=0)
plt.axhspan(0, 3.5, facecolor='red', alpha=0.2,zorder=0)
plt.xlabel(r'$L/L_\mathrm{Edd}$', fontsize=16)
plt.yticks(ticks=[1.0,2.0,3.0,4.0,5.0], labels = ['ULXP', 'SgXB', 'BeXB', 'Z-track', 'Atoll'], fontsize=16)
plt.ylim(top=5.5, bottom=0.5)
plt.xscale('log')
plt.xticks(fontsize=16)
plt.text(10,4.5,'LMXBs', fontsize=16)
plt.text(10,2.5,'HMXBs', fontsize=16)
plt.savefig('/Users/sean/Desktop/ns_eddfrac.pdf')
plt.axvline(1.0, ls='dashed', color='black')
plt.tight_layout()
plt.savefig('/Users/sean/Desktop/ns_eddfrac.pdf')