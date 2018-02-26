import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
import skimage.io as imgio
from matplotlib.offsetbox import OffsetImage, AnnotationBbox

data = fits.open('/disk/lif2/spike/detectorData/H100/H100_long_gamma_Co57_Am241_-10_gain_offset_grade0.0V.fits')[1].data

heatmap = np.reshape(data['GAIN'], [32, 32])

images = [[imgio.imread('/disk/lif2/spike/detectorData/H100/figures/pixelFits/H100_long_gamma_Co57_Am241_-10_x' + str(i) + '_y' + str(j) + '_grade0_linefit.0V.eps') for j in range(32)] for i in range(32)]

fig = plt.figure()
ax = fig.add_subplot(111)
imgplot = ax.imshow(heatmap)
cb = plt.colorbar(imgplot, label = 'Gain')
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

