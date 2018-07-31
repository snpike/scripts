import astropy.io.fits as fits
import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting

file = fits.open('/Users/sean/Desktop/SMC_X-1_final/10002013001/event_cl/nu10002013001A01_cl.evt.gz')
data = file[1].data
file.close()

countMap = np.histogram2d(data['X'], data['Y'], bins = (np.max(data['X']), np.max(data['Y'])))
edge_buffer = 10
print(np.argmax(countMap[0]))
countMap[0][:edge_buffer, :] = 0
countMap[0][:, :edge_buffer] = 0
countMap[0][-edge_buffer:, :] = 0
countMap[0][:, -edge_buffer:] = 0


x_mesh = np.array([[x for x in countMap[1][:-1]] for y in countMap[2][:-1]])
y_mesh = np.array([[y for x in countMap[1][:-1]] for y in countMap[2][:-1]])
xmax, ymax = np.unravel_index(np.argmax(countMap[0].T), np.shape(countMap[0].T))
p_init = models.Gaussian2D(amplitude = np.max(countMap[0]), x_mean = y_mesh[xmax,ymax], y_mean = x_mesh[xmax,ymax])
fit_p = fitting.LevMarLSQFitter()
p = fit_p(p_init, y_mesh, x_mesh, countMap[0].T)
print(p.x_mean)
print(p.y_mean)
CS = plt.contour(y_mesh, x_mesh, p(y_mesh, x_mesh), np.multiply(p.amplitude, np.exp([-2, -1/2])), colors = 'black')

print('done')
plt.imshow(np.log(countMap[0]))
c = plt.colorbar()
c.set_label('Counts')
plt.tight_layout()
plt.show()
plt.close()