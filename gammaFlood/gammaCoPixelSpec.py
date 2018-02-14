import matplotlib.pyplot as plt
import numpy as np
import pickle
import skimage.io as imgio
from matplotlib.offsetbox import OffsetImage, AnnotationBbox

#pixelEvents = [[[] for j in range(32)] for i in range(32)]

### Load the pickled data and do some plotting/fit the 59 keV line to a gaussian. Use this fit to convert from channels to energy.

pfile = open('/Volumes/LaCie/CdTe/gammaFlood/20170313_H100_gamma_Co57_-10C_dump.0V.pkl', 'rb')
events = pickle.load(pfile)
pfile.close()

heatmap = np.histogram2d(events['rawx'], events['rawy'], [32, 32])[0]
'''
for i in range(len(events['channel'])):
	pixelEvents[events['rawx'][i]][events['rawy'][i]].append(events['channel'][i])

del events

for x in range(32):
	for y in range(32):
		spectrum = np.histogram(pixelEvents[x][y], bins = range(3501))
		plt.figure()
		plt.plot(range(len(spectrum[0])), spectrum[0], lw = 0.5)
		plt.xlabel('Channel')
		plt.ylabel('Counts')
		plt.savefig('/Volumes/LaCie/CdTe/gammaFlood/images/pixelSpectraCo/20170313_H100_gamma_Co57_-10C_x' + str(x) + '_y' + str(y) +'_spectrum.0V.eps')
		plt.close()
'''
images = [[imgio.imread('/Volumes/LaCie/CdTe/gammaFlood/images/pixelSpectraCo/20170313_H100_gamma_Co57_-10C_x' + str(j) + '_y' + str(i) +'_spectrum.0V.eps') for j in range(32)] for i in range(32)]

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

