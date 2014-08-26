from pylab import *

from spectral import *


N = 512;
x = 2*pi/float(N);
y = 2*pi/float(N);

psi1 = rand(N,N); 
psi2 = rand(N,N); 
omega1 = zeros((N,N)); 
omega2 = zeros((N,N)); 


for i in arange(N):
    for j in arange(N):
        psi1[i,j] = sin(50*x*i)*cos(30*y*j);
        psi2[i,j] = cos(45*x*i)*sin(80*y*j);


omega1 = partialX(psi1,2) + partialY(psi1,2) + (psi2-psi1);
omega2 = partialX(psi2,2) + partialY(psi2,2) + (psi1-psi2);

diffPsi = InvPotentialVorticityTwoLayer(omega1-omega2);
sumPsi = InvLaplacian(omega1+omega2);

print 'Error in difference and sum:';
print amax(abs(diffPsi-psi1+psi2));
print amax(abs(sumPsi-psi1-psi2));
