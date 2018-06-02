import astropy.io.ascii as asciio
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

filepath = input('Please enter the directory with the leakage data: ').strip()
if filepath[-1]=='/':
	filepath = filepath[:-1]
detector = input('Please enter the detector ID: ').strip()
pos = int(input('Please enter the position of the detector: ').strip())

START = -1024 * (1 + pos)

slash = 0
i = 0
for char in filepath:
    if char == '/':
        slash = i
    i += 1

filename = filepath[slash + 1:]

Tlist  = {'-5C': [], '5C': [], '15C': [], '23C': []}
CPlist = [100, 200, 300, 400, 500, 600]
Nlist  = [300, 400, 500, 600]

outfile = open('/disk/lif2/spike/detectorData/' + detector + '/leakage.out', 'w')

for T in Tlist:
	outfile.write('T = ' + str(T) + '\n')
	# First construct the maps of ADC_0V
	CPdata = asciio.read(filepath + '/' + filename + '_' + T + '.C0V.txt')
	Ndata = asciio.read(filepath + '/' + filename + '_' + T + '.N0V.txt')
	ADC_0V_CP = np.zeros((32,32))
	ADC_0V_N = np.zeros((32,32))
	
	for i in range(1024):
		CPcol = CPdata.field('col4')[START + i]
		CProw = CPdata.field('col5')[START + i]
		Ncol = Ndata.field('col4')[START + i]
		Nrow = Ndata.field('col5')[START + i]
		ADC_0V_CP[CProw, CPcol] = CPdata.field('col6')[START + i]
		ADC_0V_N[Nrow, Ncol] = Ndata.field('col6')[START + i]
	
	for HV in CPlist:
		outfile.write('HV = ' + str(HV) + 'V\n')
		outfile.write('CP mode' + '\n')
		CPdata = asciio.read(filepath + '/' + filename + '_' + T + '.C' + str(HV) + 'V.txt')
		CPmap = np.zeros((32,32))

		for i in range(1024):
			CPcol = CPdata.field('col4')[START + i]
			CProw = CPdata.field('col5')[START + i]
			CPmap[CProw, CPcol] = (CPdata.field('col6')[START + i] - ADC_0V_CP[CProw, CPcol]) * (1.7e3)/3000

		plt.figure()
		masked = np.ma.masked_where(CPmap > 500, CPmap)
		current_cmap = mpl.cm.get_cmap()
		current_cmap.set_bad(color='gray')
		plt.imshow(masked)
		#plt.imshow(CPmap)
		c = plt.colorbar()
		c.set_label('Leakage Current (pA)')
		plt.tight_layout()
		plt.savefig('/disk/lif2/spike/detectorData/' + detector + '/figures/' + filename + '_' + T + '.C' + str(HV) + 'V.map.eps')
		plt.close()

		plt.figure()
		#plt.hist(CPmap.flatten(), bins = 50, histtype='step')
		plt.hist(masked.flatten(), bins = 50, histtype='step')
		plt.ylabel('Pixels')
		plt.xlabel('Leakage Current (pA)')
		plt.tight_layout()
		plt.savefig('/disk/lif2/spike/detectorData/' + detector + '/figures/' + filename + '_' + T + '.C' + str(HV) + 'V.hist.eps')
		plt.close()

		outfile.write('Mean leakage current: ' + str(np.mean(masked)) + '\n')
		outfile.write('leakage current standard deviation: ' + str(np.std(masked)) + '\n')
		Tlist[T].append(np.mean(masked))


		if HV in Nlist:
			outfile.write('N mode' + '\n')
			Ndata = asciio.read(filepath + '/' + filename + '_' + T  + '.N' + str(HV) + 'V.txt')
			Nmap = np.zeros((32,32))

			for i in range(1024):
				Ncol = Ndata.field('col4')[START + i]
				Nrow = Ndata.field('col5')[START + i]
				Nmap[Nrow, Ncol] = (Ndata.field('col6')[START + i] - ADC_0V_N[Nrow, Ncol]) * (1.7e3)/150

			plt.figure()
			plt.imshow(Nmap)
			c = plt.colorbar()
			c.set_label('Leakage Current (pA)')
			plt.tight_layout()
			plt.savefig('/disk/lif2/spike/detectorData/' + detector + '/figures/' + filename + '_' + T + '.N' + str(HV) + 'V.map.eps')
			plt.close()

			plt.figure()
			plt.hist(Nmap.flatten(), bins = 50, histtype='step')
			plt.ylabel('Pixels')
			plt.xlabel('Leakage Current (pA)')
			plt.tight_layout()
			plt.savefig('/disk/lif2/spike/detectorData/' + detector + '/figures/' + filename + '_' + T + '.N' + str(HV) + 'V.hist.eps')
			plt.close()

			outfile.write('Mean leakage current: ' + str(np.mean(Nmap)) + '\n')
			outfile.write('leakage current standard deviation: ' + str(np.std(Nmap)) + '\n')

outfile.close()

plt.figure()
for T in Tlist:
	plt.plot(CPlist, Tlist[T], label = r'$T=$' + T)
plt.legend()
plt.xlabel('Bias Voltage (V)')
plt.ylabel('Mean Leakage Current (pA)')
plt.tight_layout()
plt.savefig('/disk/lif2/spike/detectorData/' + detector + '/figures/' + filename + '.HV_T_plot.eps')
#plt.show()


