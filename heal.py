import matplotlib.pyplot as plt
import numpy as np
import healpy as hp
from scipy.interpolate import interp1d
#%matplotlib inline

#this function gives pdf
def gauss(x,sig):
     return np.exp(-(x**2)/(2*sig**2))

#this part of code gives the probability distribution function pdf
NSIDE = 64
print("Approximate resolution at NSIDE {} is {:.2} deg".format(NSIDE, 180*hp.nside2resol(NSIDE, arcmin=False) / 3.14))

NPIX = hp.nside2npix(NSIDE)
print NPIX
n = np.arange(0,NPIX,1)

#position of the center of the gaussian
theta0=np.pi/2
phi0=0
#gaussian parameter
sig=0.5

#computing distance and assigning probability
dist=hp.rotator.angdist(hp.pix2ang(NSIDE,n),[theta0,phi0])
m=gauss(dist,sig)
norm=np.sum(m)
m=m/norm

y0 = np.arange(NPIX)
x0 = np.zeros(NPIX)
#$\frac{12}{23}$
for i in range(0,NPIX):
    pdf=np.sum(m[0:i+1])
    # x= (x1-x0)*int + x0
    x0[i]=((pdf)*(1.) + 0) # a is a array consists of (value, error)

#print x0[-1], "surhud debug"\


f=interp1d(x0,n,kind='linear')

#sample of uniform numbers
n0=10
xnew=np.random.random_sample(n0)

ynew=f(xnew)

ynew = np.around(ynew).astype(int)
m[ynew]=np.max(m)
hp.mollview(m, title="test")
plt.show()
