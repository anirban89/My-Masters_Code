from numpy import *
from pylab import *
from scipy import *
from math import *

N = 128
omega = rand(N,N)
o = zeros((1,N))
#for i in arange(N):
#  
#       if (i<N/4):
#           omega[i,:]=0.25;
#       elif (i >= N/4 and i<3*N/4):
#	   omega[i,:]=0.5;
#       else:
#           omega[i,:]=1.5;
for i in arange(N):
	#o((i)) = i
	omega[:,i] =  0.5+(1/pi)*atan((i-N/2)/5)
#subplot(121)
imshow(omega.transpose())
#subplot(122)
#print o
#plot(o,omega[1,i])
show()

