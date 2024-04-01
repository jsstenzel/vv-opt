import numpy as np
import matplotlib.pyplot as plt
import pyds9 as ds9
from scipy import signal
from astropy.convolution import Gaussian2DKernel, convolve

naxis1 = 1024

x = np.arange(0,naxis1,1,float)
y = x[:,np.newaxis]

for i in range(3):
    for j in range(3):

        xcen = naxis1 // 4 + (i*naxis1 // 2)
        ycen = naxis1 // 4 + (j*naxis1 // 2)

#        xcen = (naxis1 // 2) + ((i-1)*291.0)
#        ycen = (naxis1 // 2) + ((j-1)*291.0)

        radius = np.sqrt((x-xcen)**2 + (y-ycen)**2)
        tophat = radius * 0
        tophat[np.where(radius < 64)] = 1.0

        if (i==0 and j==0): 
            grid = tophat
        else:
            grid += tophat

D_fib = 110.0
f_cam = 70.0
f_col = 200.0

D_geometric = D_fib * f_cam / f_col
scalefac = naxis1 / D_geometric

r_spot = 9.0 / D_geometric * 128.0 # um RMS

xcen = ycen = naxis1 // 2
radius = np.sqrt((x-xcen)**2 + (y-ycen)**2)
tophat = radius * 0
tophat[np.where(radius < 64)] = 1.0
tophat /= np.sum(tophat)

kernel = 1/np.sqrt(2*np.pi*r_spot**2) * np.exp(-1.0 * radius**2 / (2.0*r_spot**2))
kernel /= np.sum(kernel) # somehow not normalized initially

tophat_convolved = signal.fftconvolve(grid,kernel,mode='same')

## Must start DS9 from the command line, then use file->XPA->connect
## Then the communication will start to work properly.  If not, can use a method 
## from pyds9 to gram the target name if it's not DS9:ds9

targets = ds9.ds9_targets()
d = ds9.DS9(target='DS9:ds9')

master = (np.roll(tophat_convolved,4)-tophat_convolved)

grid_binned = np.zeros([naxis1/45+1,naxis1/45+1])
flex_binned = np.zeros([naxis1/45+1,naxis1/45+1])
extracted = np.zeros(naxis1/45+1)
for i in range(22):
    for j in range(22):
        grid_binned[i,j] = np.sum(tophat_convolved[i*45:(i+1)*45,j*45:(j+1)*45])/45.
        flex_binned[i,j] = np.sum(master[i*45:(i+1)*45,j*45:(j+1)*45])/45.

prof_image = np.zeros([naxis1/45+1,naxis1/45+1])
prof_image[14:20,:] = 1.0
prof_image[3:8,:] = 1.0

extracted = np.sum(prof_image*grid_binned,axis=0)
residual  = np.sum(prof_image*flex_binned,axis=0)

d.set('frame frameno 1')
d.set_np2arr(grid_binned)
d.set('zoom to fit')
d.set('scale minmax')
d.set('frame frameno 2')
d.set_np2arr(prof_image)
d.set('zoom to fit')
d.set('scale minmax')

# plt.plot(extracted)
plt.plot(residual/np.max(extracted))
plt.show()


