import numpy as np
from scipy import intergrate
#from scipy import interpolate
import matplotlib.pyplot as plt

#defination of the distribution for toy model
def func(x,a,b):#a is the mean and b is the sigma
	return np.exp(-((x-a)**2)/(2*b**2))

#creates a points for interpolation
#homogeneous(0,1) to gaussian with sigma 1 in range -1 to 1.

xnew=np.random.random_sample(100000)

#distribution parameters
mu=0.
sig=0.1
dy=sig/100

y=[]
x=[]
r0=-sig
r1=sig
norm=sp.integrate.quad(func, r0, r1, args=(mu,sig))
for i in np.arange(r0,(r1+dy),dy):
    y.append(i)
    a=sp.integrate.quad(func, r0, i, args=(mu,sig))
    # x= (x1-x0)*int + x0
    b=((a[0]/norm[0])*(1.) + 0) # a is a array consists of (value, error)
    x.append(b)
#print x

f=sp.interpolate.interp1d(x,y,kind='linear')
ynew=f(xnew)

#print xnew
#print ynew
plt.hist(xnew,label = 'x-uniform')
plt.hist(ynew,label = 'y-gaussian')
plt.legend()
plt.show()
