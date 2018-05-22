import astropy.io.ascii as asciio
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

filepath = input('Please enter the directory with the leakage data: ').strip()
if filepath[-1]=='/':
	filepath = filepath[:-1]
detector = input('Please enter the detector ID: ').strip()
pos = int(input('What is the position of the detector? ').strip())

START = -1024 * (1 + pos)

slash = 0
i = 0
for char in filepath:
    if char == '/':
        slash = i
    i += 1

filename = filepath[slash + 1:]

Tlist  = ['-5C', '5C', '15C', '23C']
CPlist = ['100V', '200V', '300V', '400V', '500V', '600V']
Nlist  = ['300V', '400V', '500V', '600V']



for T in Tlist:
	print('T = ' + str(T))
	# First construct the maps of ADC_0V
	CPdata = asciio.read(filepath + '/' + filename + '_' + T + '.C0V.txt')
	Ndata = asciio.read(filepath + '/' + filename + '_' + T + '.N0V.txt')
	ADC_0V_CP = np.zeros((32,32))
	ADC_0V_N = np.zeros((32,32))
	
	for i in range(1024):
		CPx = CPdata.field('col4')[START + i]
		CPy = CPdata.field('col5')[START + i]
		Nx = Ndata.field('col4')[START + i]
		Ny = Ndata.field('col5')[START + i]
		ADC_0V_CP[CPx, CPy] = CPdata.field('col6')[START + i]
		ADC_0V_N[Nx, Ny] = Ndata.field('col6')[START + i]
	
	for HV in CPlist:
		print('HV = ' + str(HV))
		print('CP mode')
		CPdata = asciio.read(filepath + '/' + filename + '_' + T + '.C' + HV + '.txt')
		CPmap = np.zeros((32,32))

		for i in range(1024):
			CPx = CPdata.field('col4')[START + i]
			CPy = CPdata.field('col5')[START + i]
			CPmap[CPx, CPy] = (CPdata.field('col6')[START + i] - ADC_0V_CP[CPx, CPy]) * (1.7e3)/3000

		plt.figure()
		plt.imshow(CPmap)
		c = plt.colorbar()
		c.set_label('Leakage Current (pA)')
		plt.tight_layout()
		plt.savefig('/disk/lif2/spike/detectorData/' + detector + '/figures/' + filename + '_' + T + '.C' + HV + '.map.eps')
		plt.close()

		plt.figure()
		plt.hist(CPmap.flatten(), bins = 50, histtype='step')
		plt.ylabel('Pixels')
		plt.xlabel('Leakage Current (pA)')
		plt.tight_layout()
		plt.savefig('/disk/lif2/spike/detectorData/' + detector + '/figures/' + filename + '_' + T + '.C' + HV + '.hist.eps')
		plt.close()

		print('Mean leakage current: ' + str(np.mean(CPmap)))
		print('leakage current standard deviation: ' + str(np.std(CPmap)))


		if HV in Nlist:
			print('N mode')
			Ndata = asciio.read(filepath + '/' + filename + '_' + T  + '.N' + HV + '.txt')
			Nmap = np.zeros((32,32))

			for i in range(1024):
				Nx = Ndata.field('col4')[START + i]
				Ny = Ndata.field('col5')[START + i]
				Nmap[Nx, Ny] = (Ndata.field('col6')[START + i] - ADC_0V_N[Nx, Ny]) * (1.7e3)/150

			plt.figure()
			plt.imshow(Nmap)
			c = plt.colorbar()
			c.set_label('Leakage Current (pA)')
			plt.tight_layout()
			plt.savefig('/disk/lif2/spike/detectorData/' + detector + '/figures/' + filename + '_' + T + '.N' + HV + '.map.eps')
			plt.close()

			plt.figure()
			plt.hist(Nmap.flatten(), bins = 50, histtype='step')
			plt.ylabel('Pixels')
			plt.xlabel('Leakage Current (pA)')
			plt.tight_layout()
			plt.savefig('/disk/lif2/spike/detectorData/' + detector + '/figures/' + filename + '_' + T + '.N' + HV + '.hist.eps')
			plt.close()

			print('Mean leakage current: ' + str(np.mean(Nmap)))
			print('leakage current standard deviation: ' + str(np.std(Nmap)))