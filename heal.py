# This code is used to produce a GW source from  a pdf.
# for the toy model we are using a gaussian centered around DM particles with a constant sigma.
# This code takes two inputs file with positions of DM particles and a sigma
import matplotlib.pyplot as plt
import numpy as np
import healpy as hp
from scipy.interpolate import interp1d

# Please note that here the positions use theta(angle from z),phi(angle from x).
# and all are in radians

# For info related to healpix in general please look at this
# https://arxiv.org/abs/astro-ph/0409513

# for more info related to healpy functions please open the link below
# https://healpy.readthedocs.io/en/latest/healpy_pix.html#conversion-from-to-sky-coordinates

#This function gives pdf
def gauss(x,sig):
     return np.exp(-(x**2)/(2*sig**2))

#this part of code gives the probability distribution function pdf
# NSIDE is the parameter for healpix structure
NSIDE = 32
print("Approximate resolution at NSIDE {} is {:.2} deg".format(NSIDE, 180*hp.nside2resol(NSIDE, arcmin=False) / 3.14))

# hp.nside2npix() computes total number of pixels = 12*NSIDE**2 
NPIX = hp.nside2npix(NSIDE)
n = np.arange(0,NPIX,1)

# Here we are loading data 
# Data(s.no,D_c,D_l,r.a,dec) - data structure
# s.no - serial no
# D_c - comoving distance
# D_l - luminosity distance
# r_a - right accertion
# dec - declination
data =np.loadtxt('/Users/divyarana/Documents/surhud/data/angular_DM.dat',dtype=float,unpack=True)
nmax=len(data[0][:])

# opening new file to write
lb=open('new_positions.txt','w')

# In this part of the code we will produce a new position provided old ones from the data
for i in range (0,nmax):
	#position of the center of the gaussian
	theta0=np.pi/2 - data[4][i] #declination
	phi0=data[3][i] #right acertion

	#gaussian parameter sig should be choosen according to the Theta max in the cross angular correlation 
	sig=0.5 #(data[2][i]*np.pi)/(np.sqrt(np.pi)*180)  # pi*sig^2/dl^2 = 1 str

	#computing distance and assigning probability
	# pix2ang : nside,ipix -> theta[rad](0,pi),phi[rad](0,2pi) 
	# hp.rotator.angdist(dir1,dir2) -> Returns the angular distance between dir1 and dir2.
	dist=hp.rotator.angdist(hp.pix2ang(NSIDE,n),[theta0,phi0])
	m=gauss(dist,sig)
	norm=np.sum(m)
	m=m/norm
	x0 = np.zeros(NPIX)
	for i1 in range(0,NPIX):
    		pdf=np.sum(m[0:i1+1])
    		# x= (x1-x0)*int + x0
    		x0[i1]=((pdf)*(1.) + 0) # a is a array consists of (value, error)

	f=interp1d(x0,n,kind='linear')

	#sample of uniform numbers
	n0=1
	xnew=np.random.random_sample(n0)
	ynew=f(xnew)
	ynew = np.around(ynew).astype(int)
	new_pos = hp.pix2ang(NSIDE,ynew[0])
	print >> lb,'',data[0][i],'\t',data[1][i],'\t',data[2][i],'\t',data[3][i],'\t',data[4][i],'\t',new_pos[1],'\t',np.pi/2 - new_pos[0]
	
	print'particle number = %d done ' % (i) # This will print the running status of the code
lb.close()

print 'New position written in new_position.txt'

m[ynew[0]]=np.max(m)
hp.mollview(m, title="test")
plt.show()

