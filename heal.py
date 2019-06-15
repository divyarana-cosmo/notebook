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

#data(s.no,D_c,D_l,r.a,dec)
data =np.loadtxt('/Users/divyarana/Documents/surhud/data/angular_DM.dat',dtype=float,unpack=True)

nmax=len(data[0][:])

lb=open('new_positions.txt','w')
#looping over all sources
for i in range (0,1):

	#position of the center of the gaussian
	theta0=np.pi/2 - data[4][i] #declination
	phi0=data[3][i] #right acertion

	#gaussian parameter
	sig=0.5#(data[2][i]*np.pi)/(np.sqrt(np.pi)*180)  # pi*sig^2/dl^2 = 1 str

	#computing distance and assigning probability
	dist=hp.rotator.angdist(hp.pix2ang(NSIDE,n),[theta0,phi0])
	m=gauss(dist,sig)
	norm=np.sum(m)
	m=m/norm

	x0 = np.zeros(NPIX)
	#$\frac{12}{23}$
	for i1 in range(0,NPIX):
    		pdf=np.sum(m[0:i1+1])
    		# x= (x1-x0)*int + x0
    		x0[i1]=((pdf)*(1.) + 0) # a is a array consists of (value, error)

	f=interp1d(x0,n,kind='linear')

	#sample of uniform numbers
	n0=1000
	xnew=np.random.random_sample(n0)
	ynew=f(xnew)
	ynew = np.around(ynew).astype(int)
	new_pos = hp.pix2ang(NSIDE,ynew)

	print >> lb,'',data[0][i],'\t',data[1][i],'\t',data[2][i],'\t',data[3][i],'\t',data[4][i],'\t',new_pos[1],'\t',np.pi/2 - new_pos[0]
	print i
lb.close()

print 'writing done'

m[ynew]=np.max(m)
hp.mollview(m, title="test")
plt.show()


'''
ynew = np.around(ynew).astype(int)
m[ynew]=np.max(m)
hp.mollview(m, title="test")
plt.show()
'''
