from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import pickle
import skimage.io as imgio
from matplotlib.offsetbox import OffsetImage, AnnotationBbox

pixelEvents = [[[] for j in range(32)] for i in range(32)]

file = fits.open('/Volumes/LaCie/CdTe/gammaFlood/20170313_H100_gamma_Co57_-10C.0V.fits')
data = file[1].data
file.close()
START = 0
i = 0

# Skip the beginning
while data.field('TEMP')[i] < -50:
    i += 1
START = i

'''
while i < len(data.field('PH')):
    if not(data.field('STIM')[i]):
        rawx = data.field('RAWX')[i]
        rawy = data.field('RAWY')[i]
        temp = data.field('PH_COM')[i].reshape(3,3)
        if np.sum(temp) > 0 and (data.field('GRADE')[i]==0):
            mask = (temp > 0).astype(int)
            channel = np.sum(np.multiply(mask, temp))
            if (not np.isnan(channel)):
                pixelEvents[rawx][rawy].append(channel)
    i += 1

for x in range(32):
	for y in range(32):
		spectrum = np.histogram(pixelEvents[x][y], bins = range(12500))
		plt.figure()
		plt.plot(range(len(spectrum[0])), spectrum[0], lw = 0.5)
		plt.xlabel('Channel')
		plt.ylabel('Counts')
		plt.savefig('/Volumes/LaCie/CdTe/gammaFlood/images/pixelSpectraCoUncor/20170313_H100_gamma_Co57_-10C_x' + str(x) + '_y' + str(y) +'_spectrum.0V.eps')
		plt.close()

'''
delList = []
i = START
while i < len(data.field('PH')):
    temp = data.field('PH_COM')[i].reshape(3,3)
    if data.field('STIM')[i] or (np.sum(temp) <= 0 and (data.field('GRADE')[i]!=0)):
        delList.append(i-START)
    i += 1

heatmap = np.histogram2d(np.delete(data.field('RAWX')[START:], delList), np.delete(data.field('RAWY')[START:], delList), [32, 32])[0]

images = [[imgio.imread('/Volumes/LaCie/CdTe/gammaFlood/images/pixelSpectraCoUncor/20170313_H100_gamma_Co57_-10C_x' + str(j) + '_y' + str(i) +'_spectrum.0V.eps') for j in range(32)] for i in range(32)]

fig = plt.figure()
ax = fig.add_subplot(111)
imgplot = ax.imshow(heatmap)
cb = plt.colorbar(imgplot, label = 'Counts')
cb.ax.zorder = -1

# create the annotations box
im = OffsetImage(images[0][0], zoom = 0.7)
xybox=(150., 150.)
ab = AnnotationBbox(im, (0,0), xybox=xybox, xycoords='data',
        boxcoords="offset points",  pad=0.1,  arrowprops=dict(arrowstyle="->"))
# add it to the axes and make it invisible
ax.add_artist(ab)
ab.set_visible(False)

def hover(event):
    # if the mouse is over the scatter points
    if imgplot.contains(event)[0]:
        # get the figure size
        #w,h = fig.get_size_inches()*fig.dpi
        ws = (event.xdata > 15.5)*-1 + (event.xdata <= 15.5) 
        hs = (event.ydata > 15.5) + (event.ydata <= 15.5)*-1
        # if event occurs in the top or right quadrant of the figure,
        # change the annotation box position relative to mouse.
        ab.xybox = (xybox[0]*ws, xybox[1]*hs)
        # make annotation box visible
        ab.set_visible(True)
        # place it at the position of the hovered scatter point
        ab.xy =(event.xdata, event.ydata)
        # set the image corresponding to that point
        im.set_data(images[int(round(event.xdata))][int(round(event.ydata))])
    else:
        #if the mouse is not over a scatter point
        ab.set_visible(False)
    fig.canvas.draw_idle()

# add callback for mouse moves
fig.canvas.mpl_connect('motion_notify_event', hover)      
plt.show()

