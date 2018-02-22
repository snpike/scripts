from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np


formats = {'TIME': 'D', 'CHANNEL': 'D', 'GRADE': 'K', 'STIM': 'K', 'PH_COM': '9D', 'SOURCE': '5A'}

for source in [['Co57','/disk/lif2/spike/detectorData/H100/longGammaFlood/20170908_H100_long_gamma_Co57_-10.0V.fits'], ['Am241', '/disk/lif2/spike/detectorData/H100/longGammaFlood/20170913_H100_long_gamma_Am241_-10.0V.fits']]:

    file = fits.open(source[1], memmap=True)
    data = file[1].data

    sortedIndices = np.argsort(data, order=('RAWX', 'RAWY'))

    i = 0
    for x in range(32):
        print(str(x))
        for y in range(32):
            collist = []
            newdata = {'TIME': [], 'CHANNEL': [], 'GRADE': [], 'STIM': [], 'PH_COM': [], 'SOURCE': []}
            while (data.field('RAWY')[sortedIndices[i]]==y):
                if data.field('TEMP')[sortedIndices[i]] > -50:
                    temp = data.field('PH_COM')[sortedIndices[i]].reshape(3,3)
                    if np.sum(temp) > 0:
                        mask = (temp > 0).astype(int)
                        channel = np.sum(np.multiply(mask, temp))
                        if (not np.isnan(channel)):
                            newdata['TIME'].append(data.field('TIME')[sortedIndices[i]])
                            newdata['CHANNEL'].append(channel)
                            newdata['GRADE'].append(data.field('GRADE')[sortedIndices[i]])
                            newdata['STIM'].append(data.field('STIM')[sortedIndices[i]])
                            newdata['PH_COM'].append(data.field('PH_COM')[sortedIndices[i]])
                            newdata['SOURCE'].append(source[0])
                i += 1
                if not (i < len(sortedIndices)):
                    break
            for key in newdata:
                collist.append(fits.Column(name=key, format=formats[key], array=newdata[key]))
            fits.BinTableHDU.from_columns(collist).writeto('/disk/lif2/spike/detectorData/H100/longGammaFlood/pixelData/H100_long_gamma_' + source[0] + '_-10_0V_x' + str(x) + '_y' + str(y) + '.fits')

print('done')


'''
    for i in sortedIndices:
        temp = data.field('PH_COM')[i].reshape(3,3)
        if np.sum(temp) > 0:
            mask = (temp > 0).astype(int)
            channel = np.sum(np.multiply(mask, temp))
            if (not np.isnan(channel)):
                tmpfile=fits.open(files[data.field('RAWX')[i], data.field('RAWY')[i]], memmap=True, mode='update')
                tmpfile[1].data['TIME'] = np.append(tmpfile[1].data['TIME'], data.field('TIME')[i])
                tmpfile[1].data['CHANNEL'] = np.append(tmpfile[1].data['CHANNEL'], channel)
                tmpfile[1].data['GRADE'] = np.append(tmpfile[1].data['GRADE'], data.field('GRADE')[i])
                tmpfile[1].data['STIM'] = np.append(tmpfile[1].data['STIM'], data.field('STIM')[i])
                tmpfile[1].data['PH_COM'] = np.append(tmpfile[1].data['PH_COM'], data.field('PH_COM')[i])
                tmpfile[1].data['SOURCE'] = np.append(tmpfile[1].data['SOURCE'], source[0])
                tmpfile.flush()
                tmpfile.close()
        i += 1

# Apply the first-level gain and offset corrections
while i < len(data.field('PH_COM')):
	rawx = data.field('RAWX')[i]
	rawy = data.field('RAWY')[i]
	temp = data.field('PH_COM')[i].reshape(3,3)
	if np.sum(temp) > 0:
		mask = (temp > 0).astype(int)
		channel = np.sum(np.add(np.multiply(np.multiply(mask, temp), gain[rawx:rawx + 3, rawy:rawy + 3]), np.multiply(mask, offset[rawx:rawx + 3, rawy:rawy + 3])))
		if (not np.isnan(channel)):
			events['RAWX'].append(rawx)
			events['RAWY'].append(rawy)
			events['TIME'].append(data.field('TIME')[i])
			events['CHANNEL'].append(channel)
			events['GRADE'].append(data.field('GRADE')[i])
			events['STIM'].append(data.field('STIM')[i])
			events['PH_COM'].append(temp)
	i += 1


spectrum = np.histogram(events['CHANNEL'], bins=int(np.ceil(np.max(events['CHANNEL']))))
plt.figure()
plt.plot(range(len(spectrum[0])), spectrum[0])
plt.show()
plt.close()'''

