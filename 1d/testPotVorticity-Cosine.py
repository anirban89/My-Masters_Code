from spectral import *;
from pylab import *;

N = 2048;
ly = 8*pi;
lx = 4*ly;
t = linspace(0,ly,N/2+1);

tx = linspace(0,lx,N+1);
tx = tx[0:len(tx)-1];

test = zeros((N,N/2+1));
print 'sin', shape(sin(t)), 'test', shape(test);

for i in arange(N):
    test[i,:] = sin(tx[i])*cos(t);


out = InvPotentialVorticityCosine(test);

print shape(out);

new = partialYCosine(out,order=2) + \
    partialX(out,order=2) - out;


print 'max value in test: ', amax(test);
print 'max value in out: ', amax(out);
print 'max value in new: ', amax(new);
imshow(test);
colorbar();
figure();
imshow(new);
colorbar();
figure();
imshow(new-test);
colorbar();
show();


