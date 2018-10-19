from astropy.io import fits
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.modeling import models, fitting
import os.path
import pickle
import seaborn as sns
from matplotlib.gridspec import GridSpec
import scipy.stats as stats

sns.set_context('talk')
sns.set_style("ticks")
sns.set_palette("colorblind")

# Line data from https://www-nds.iaea.org/relnsd/vcharthtml/VChartHTML.html

# sourcelist = {'Am241': [[59.541], r'${}^{241}{\rm Am}$', input('Please enter the filepath to the Am241 flood data: ').strip()], \
#             'Co57': [[14.41, 122.06], r'${}^{57}{\rm Co}$', input('Please enter the filepath to the Co57 flood data: ').strip()], \
#             'Eu155': [[86.548], r'${}^{155}{\rm Eu}$', input('Please enter the filepath to the Eu155 flood data: ').strip()]}

sourcelist = ['Am241', 'Co57', 'Eu155']

line_dict = {'Am241': [59.541], 'Co57': [14.41, 122.06], 'Eu155': [86.548]}
latex_dict = {'Am241': r'${}^{241}{\rm Am}$', 'Co57': r'${}^{57}{\rm Co}$', 'Eu155': r'${}^{155}{\rm Eu}$'}
filepath_dict = {'Am241': input('Please enter the filepath to the Am241 flood data: ').strip(), \
                'Co57': input('Please enter the filepath to the Co57 flood data: ').strip(), \
                'Eu155': input('Please enter the filepath to the Eu155 flood data: ').strip()}

for source in sourcelist:
    if not os.path.exists(sourcelist[source][2]):
        print('Could not find the file entered for ' + source + '. The script will now exit.')
        exit()

gainpath = input('Please enter the filepath to the gain data: ').strip()

detector = input('Please enter the detector ID: ').strip()

temp_response = input('Should I analyze the full detector? (y/n) ').strip()
while temp_response not in ['y','yes', 'n', 'no']:
    temp_response = input('Should I analyze the full detector? (y/n) ').strip()

full_det = (temp_response in ['y','yes'])
if full_det:
    region = [[0,32][0,32]]
else:
    region = [[int(input('Low column? ').strip()), int(input('Low row? ').strip())],[int(input('High column? ').strip()), int(input('High row? ').strip())]]

gainBool = os.path.exists(gainpath)

maxchannel = 10000
bins = 10000

spectra = {}

buff_gain = np.zeros((34, 34))

if gainBool:
    buff_gain[1:33][1:33] = numpy.loadtxt(gainpath)

else:
    # Collect energy-channel pairs at each spectral line which we will later fit to a polynomial to find the gain
    gain_points = [[[] for x in range(32)] for y in range(32)]

    gain = np.zeros((32, 32))

    for source in sourcelist:
        lines = line_dict[source]
        latex_label = latex_dict[source]
        filepath = filepath_dict[source]

        slash = 0
        i = 0
        for char in filepath:
            if char == '/':
                slash = i
            i += 1

        filename = filepath[slash + 1:]

        file = fits.open(filepath)
        data = file[1].data
        file.close()

        T_mask = data['TEMP'] > -20

        START = np.argmax(T_mask)
        END = len(T_mask) - np.argmax(T_mask[::-1])

        STIMmask = np.array(data['STIM'][START:END])==0

        for row in np.arange(region[0][1], region[1][1] + 1):

            row_mask = data['RAWY'] == row
            
            for col in np.arange(region[0][0], region[1][0] + 1):
                plt.figure()

                col_mask = data['RAWX'] == col
                channel = data.field('PH')[START:END][np.nonzero(np.multiply(np.multiply((col_mask), (row_mask)), data.field('GRADE')[START:END] == 0))]
            
                if len(channel):

                    spectrum = plt.hist(channel, bins=bins, range = (0, maxchannel))
                    
                    for line in lines:
                    
                        centroid = np.argmax(spectrum[0][int(line/0.013)-500:int(line/0.013)+500]) + int(line/0.013)-500
                        fit_channels = np.arange(centroid-70, centroid + 150)
                        g_init = models.Gaussian1D(amplitude=spectrum[0][centroid], mean=centroid, stddev = 75)
                        fit_g = fitting.LevMarLSQFitter()
                        g = fit_g(g_init, fit_channels, spectrum[0][fit_channels])
                
                        if fit_g.fit_info['param_cov'] is not None:
                            
                            plt.plot(fit_channels, g(fit_channels), label = 'Gaussian fit (' + str(round(line, 0)) + ' keV)')
                            plt.ylabel('Counts')
                            plt.xlabel('Channel')
                            plt.legend()
                            plt.savefig('/users/spike/det_figs/' + detector + '/pixels/' + filename[:-5] + '_x' + str(x) + '_y' + str(y) + '_spec.pdf')

                            # If the gain is not near  0.013 then the spectrum was probably not good enough to get a real gain value. Skip it
                            if 0.01 < (line/g.mean) < 0.016:
                                # Put the line energy (in keV) and the mean in channels in the gain data points
                                gain_points[row][col].append([g.mean, line])

    for row in range(32):
        for col in range(32):
            if len(gain_points[row][col]):
                points = np.array(gain_points[row][col]).T
                gain[row][col] = stats.linregress(points[0], y = points[1])[0]

    # interpolate gain for pixels where fit was unsuccessful
    newgain = np.zeros((34,34))
    newgain[1:33, 1:33] = gain
    empty = np.transpose(np.nonzero(gain == 0.0))
    #print(empty)
    for x in empty:
        if region[0][1] < x[0] < region[1][1]:
            if region[0][0] < x[1] < region[1][0]:
                temp = newgain[x[0]:x[0]+3, x[1]:x[1]+3]
                if np.count_nonzero(temp):
                    gain[x[0], x[1]] = np.sum(temp)/np.count_nonzero(temp)

    # One more pass for good measure
    newgain = np.zeros((34,34))
    newgain[1:33, 1:33] = gain
    empty = np.transpose(np.nonzero(gain == 0.0))
    #print(empty)
    for x in empty:
        if region[0][1] < x[0] < region[1][1]:
            if region[0][0] < x[1] < region[1][0]:
                temp = newgain[x[0]:x[0]+3, x[1]:x[1]+3]
                if np.count_nonzero(temp):
                    gain[x[0], x[1]] = np.sum(temp)/np.count_nonzero(temp)

    np.savetxt('/disk/lif2/spike/detectorData/' + detector + '/fullgain_region_low_x' + str(region[0][0]) + '_y' + str(region[0][1]) + 'high_x' + str(region[1][0]) + '_y' + str(region[1][1]) + '.npy', gain)
    buff_gain[1:33][1:33] = gain


for source in sourcelist:

    lines = line_dict[source]
    latex_label = latex_dict[source]
    filepath = filepath_dict[source]

    slash = 0
    i = 0
    for char in filepath:
        if char == '/':
            slash = i
        i += 1

    filename = filepath[slash + 1:]

    file = fits.open(filepath)
    data = file[1].data
    file.close()

    T_mask = data['TEMP'] > -20

    START = np.argmax(T_mask)
    END = len(T_mask) - np.argmax(T_mask[::-1])

    STIMmask = np.array(data['STIM'][START:END])==0

    countMap = [[np.sum(np.multiply(STIMmask, np.multiply(np.array(data['RAWX'][START:END])==col, np.array(data['RAWY'][START:END])==row))) for col in np.arange(region[0][0], region[1][0] + 1)] for row in np.arange(region[0][1], region[1][1] + 1)]

    count_array = np.array(countMap)

    # Plot some info about counts
    plt.figure()
    masked = np.ma.masked_values(countMap, 0.0)
    current_cmap = mpl.cm.get_cmap()
    current_cmap.set_bad(color='gray')
    plt.imshow(masked)
    c = plt.colorbar()
    c.set_label('Counts')
    plt.xticks(np.arange(0,region[1][0]-region[0][0]+1), [str(int(x)) for x in np.arange(region[0][0], region[1][0] + 1)])
    plt.yticks(np.arange(0,region[1][1]-region[0][1]+1), [str(int(x)) for x in np.arange(region[0][1], region[1][1] + 1)])
    plt.tight_layout()
    plt.savefig('/users/spike/det_figs/' + detector + '/' + filename[:-5] + '_floodmap.pdf')
    plt.close()

    plt.figure()
    plt.hist(np.array(countMap).flatten(), bins = 100, range = (0, np.max(countMap) + 1), histtype = 'stepfilled')
    plt.ylabel('Pixels')
    plt.xlabel('Counts')
    plt.tight_layout()
    plt.savefig('/users/spike/det_figs/' + detector + '/' + filename[:-5] + '_gammahist.pdf')
    plt.close()

    # If there's gain data then correct the spectrum
    energyList = []


    for row in np.arange(region[0][1], region[1][1] + 1):
        row_mask = data['RAWY'] == row
        for col in np.arange(region[0][0], region[1][0] + 1):
            col_mask = data['RAWX'] == col
            # Getting indices ('inds') and PH_COM values ('pulses') of 
            # all events at current pixel.
            inds = np.nonzero(np.multiply(np.multiply(row_mask, col_mask), T_mask))
            pulses = data.field('PH_COM')[inds]
            # The gain for the 3x3 grid around this pixel
            gain_grid = buff_gain[row:row + 3, col:col + 3]
            # iterating through the PH_COM values for this pixel
            for pulse in pulses:
                # Append the sum of positive energies in the 
                # pulse grid to 'energies'
                pulse_grid = pulse.reshape(3, 3)
                mask = (pulse_grid > 0).astype(int)
                energyList.append(np.sum(np.multiply(np.multiply(mask, pulse_grid), gain_grid)))


    bins = 10000
    spectra[source] = np.histogram(energyList, bins = bins, range= (0.01, 120))
    plt.figure()
    plt.plot(spectrum[1][:-1], spectrum[0], label = latex_label)

    for line in lines:
        centroid = np.argmax(spectrum[0][int(line/0.013)-500:int(line/0.013)+500]) + int(line/0.013)-500
        fit_channels = np.arange(centroid-70, centroid + 150)
        g_init = models.Gaussian1D(amplitude=spectrum[0][centroid], mean=centroid, stddev = 75)
        fit_g = fitting.LevMarLSQFitter()
        g = fit_g(g_init, fit_channels, spectrum[0][fit_channels])
        sigma_err = np.diag(fit_g.fit_info['param_cov'])[2]
        fwhm_err = 2*np.sqrt(2*np.log(2))*sigma_err
        mean_err = np.diag(fit_g.fit_info['param_cov'])[1]
        frac_err = np.sqrt(np.square(fwhm_err) + np.square(g.fwhm*mean_err/g.mean))/g.mean
        plt.plot(spectrum[1][fit_channels], g(fit_channels), label = 'Gaussian fit')
        print('FWHM ' + str(g.fwhm) + '+/-' + str(fwhm_err) + ' at ' str(round(line, 0)) + ' keV')

    plt.xlabel('Energy (keV)')
    plt.ylabel('Counts')
    plt.legend()
    plt.tight_layout()
    plt.savefig('/users/spike/det_figs/' + detector + '/' + filename[:-5] + '_spec_gain.pdf')
    plt.close()

plt.figure()

for source in sourcelist:
    lines = line_dict[source]
    latex_label = latex_dict[source]
    spectrum = spectra[source]

    plt.plot(spectrum[1][:-1], spectrum[0], label = latex_label)

    for line in lines:
        centroid = np.argmax(spectrum[0][int(line/0.013)-500:int(line/0.013)+500]) + int(line/0.013)-500
        fit_channels = np.arange(centroid-70, centroid + 150)
        g_init = models.Gaussian1D(amplitude=spectrum[0][centroid], mean=centroid, stddev = 75)
        fit_g = fitting.LevMarLSQFitter()
        g = fit_g(g_init, fit_channels, spectrum[0][fit_channels])
        sigma_err = np.diag(fit_g.fit_info['param_cov'])[2]
        fwhm_err = 2*np.sqrt(2*np.log(2))*sigma_err
        mean_err = np.diag(fit_g.fit_info['param_cov'])[1]
        frac_err = np.sqrt(np.square(fwhm_err) + np.square(g.fwhm*mean_err/g.mean))/g.mean
        plt.plot(spectrum[1][fit_channels], g(fit_channels), label = 'Gaussian fit')

plt.xlabel('Energy (keV)')
plt.ylabel('Counts')
plt.legend()
plt.tight_layout()
plt.savefig('/users/spike/det_figs/' + detector + '/fullspec_region_low_x' + str(region[0][0]) + '_y' + str(region[0][1]) + 'high_x' + str(region[1][0]) + '_y' + str(region[1][1]) + '.pdf')
plt.close()


