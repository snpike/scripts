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

spectra = {}

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

    pix_on = np.nonzero(countMap)
    count_array = np.array(countMap)

    # Plot some info about counts
    plt.figure()
    masked = np.ma.masked_values(countMap, 0.0)
    current_cmap = mpl.cm.get_cmap()
    current_cmap.set_bad(color='gray')
    plt.imshow(masked)
    c = plt.colorbar()
    c.set_label('Counts')
    plt.xticks(np.arange(0,region[1][0]-region[0][0]+1), np.arange(region[0][0], region[1][0] + 1))
    plt.yticks(np.arange(0,region[1][1]-region[0][1]+1), np.arange(region[0][1], region[1][1] + 1))
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
    resetList = []
    priorList = []
    if gainBool:
        gain = np.zeros((34, 34))
        gain[1:33,1:33] = pickle.load(open(gainpath, 'rb'))
        # for event in data[START:END]:
        #     row = event['RAWY']
        #     col = event['RAWX']
        #     if [row, col] not in badpix:
        #         temp = event['PH_COM'].reshape(3,3)
        #         mask = (temp > 0).astype(int)
        #         energyList.append(np.sum(np.multiply(np.multiply(mask, temp), gain[row:row + 3, col:col + 3])))

        for row in np.arange(region[0][1], region[1][1] + 1):
            row_mask = data['RAWY'] == row
            for col in np.arange(region[0][0], region[1][0] + 1):
                if [row, col] not in badpix:
                    col_mask = data['RAWX'] == col
                    # Getting indices ('inds') and PH_COM values ('pulses') of 
                    # all events at current pixel.
                    inds = np.nonzero(np.multiply(np.multiply(row_mask, col_mask), T_mask))
                    pulses = data.field('PH_COM')[inds]
                    # The gain for the 3x3 grid around this pixel
                    gain_grid = gain[row:row + 3, col:col + 3]
                    priorList = priorList + list(data.field('PRIOR')[inds])
                    resetList = resetList + list(data.field('RESET')[inds])
                    # iterating through the PH_COM values for this pixel
                    for pulse in pulses:
                        # Append the sum of positive energies in the 
                        # pulse grid to 'energies'
                        pulse_grid = pulse.reshape(3, 3)
                        mask = (pulse_grid > 0).astype(int)
                        energyList.append(np.sum(np.multiply(
                            np.multiply(mask, pulse_grid), gain_grid)))


        bins = 10000
        spectra[source] = np.histogram(energyList, bins = bins, range= (0.01, 120))

        for line in lines:
            centroid = np.argmax(spectrum[0][int(line/0.013)-500:int(line/0.013)+500]) + int(line/0.013)-500
            fit_channels = np.arange(centroid-70, centroid + 150)
            g_init = models.Gaussian1D(amplitude=spectrum[0][centroid], mean=centroid, stddev = 75)
            fit_g = fitting.LevMarLSQFitter()
            g = fit_g(g_init, fit_channels, spectrum[0][fit_channels])
            # print(np.diag(fit_g.fit_info['param_cov']))
            sigma_err = np.diag(fit_g.fit_info['param_cov'])[2]
            fwhm_err = 2*np.sqrt(2*np.log(2))*sigma_err
            mean_err = np.diag(fit_g.fit_info['param_cov'])[1]
            frac_err = np.sqrt(np.square(fwhm_err) + np.square(g.fwhm*mean_err/g.mean))/g.mean
            # print(g.fwhm/g.mean)
            # print(frac_err)
            # print(Am_line * g.fwhm/g.mean)
            # print(frac_err * Am_line)
            plt.text(70, spectrum[0][centroid]*3/5, r'$\mathrm{FWHM}=$' + str(int(round(Am_line * 1000 * g.fwhm/g.mean, 0))) + r'$\pm$' + str(int(round(frac_err * Am_line*1000))) + ' eV', fontsize=13)

            plt.plot(spectrum[1][:-1], spectrum[0], label = r'${}^{241}{\rm Am}$')
            plt.plot(spectrum[1][fit_channels], g(fit_channels), label = 'Gaussian fit')
            plt.xlabel('Energy (keV)')
            plt.ylabel('Counts')
            plt.legend()

            plt.tight_layout()
            plt.savefig('/users/spike/det_figs/' + detector + '/' + filename[:-5] + '_gammaspec_gain.pdf')
            plt.close()

    else:
        bins = np.arange(1,maxchannel)
        spectrum = np.histogram(data['PH'][START:END], bins = bins, range= (0, maxchannel))

        centroid = np.argmax(spectrum[0][3000:]) + 3000
        fit_channels = np.arange(centroid-100, centroid + 200)
        g_init = models.Gaussian1D(amplitude=spectrum[0][centroid], mean=centroid, stddev = 75)
        fit_g = fitting.LevMarLSQFitter()
        g = fit_g(g_init, fit_channels, spectrum[0][fit_channels])
        print(np.diag(fit_g.fit_info['param_cov']))
        sigma_err = np.diag(fit_g.fit_info['param_cov'])[2]
        fwhm_err = 2*np.sqrt(2*np.log(2))*sigma_err
        mean_err = np.diag(fit_g.fit_info['param_cov'])[1]
        frac_err = np.sqrt(np.square(fwhm_err) + np.square(g.fwhm*mean_err/g.mean))/g.mean
        print(g.fwhm/g.mean)
        print(frac_err)
        print(Am_line * g.fwhm/g.mean)
        print(frac_err * Am_line)
        #plt.text(maxchannel*3/5, spectrum[0][centroid]*3/5, r'$\mathrm{FWHM}=$' + str(int(g.fwhm)) + r'$\pm$' + str(int(2*np.sqrt(2*np.log(2))*sigma_err)), fontsize=16)
        #plt.text(maxchannel*3/5, spectrum[0][centroid]*3/5, r'$\mathrm{\frac{FWHM}{\mu}}=$' + str(int(round(100*g.fwhm/g.mean, 0))) + '%', fontsize=14)
        plt.text(maxchannel*3/5, spectrum[0][centroid]*3/5, r'$\mathrm{FWHM}=$' + str(int(round(Am_line * 1000 * g.fwhm/g.mean, 0))) + r'$\pm$' + str(int(round(frac_err * Am_line*1000))) + ' eV', fontsize=13)

        plt.plot(spectrum[1][:-1], spectrum[0], label = r'${}^{241}{\rm Am}$')
        plt.plot(fit_channels, g(fit_channels), label = 'Gaussian fit')
        plt.xlabel('Channel')
        plt.ylabel('Counts')
        plt.legend()

        #plt.title(detector + ' ' + source + ' Spectrum ' + '(' + etc + ')')
        plt.tight_layout()
        #plt.savefig('/disk/lif2/spike/detectorData/' + detector + '/figures/' + filename[:-4] + 'gammaspec.pdf')
        plt.savefig('/users/spike/det_figs/' + detector + '/' + filename[:-5] + '_gammaspec.pdf')
        #plt.show()
        plt.close()