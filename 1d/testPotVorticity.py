from spectral import *;
from pylab import *;

N = 1280;
t = linspace(0,8*pi,N+1);
t = t[0:len(t)-1];

test = zeros((N,N));
print 'sin', shape(sin(t)), 'test', shape(test);

for i in arange(N):
    test[i,:] = sin(t[i])*cos(t);


out = InvPotentialVorticity(test,length=2*pi);

new = partialY(out,order=2,length=2*pi) + \
    partialX(out,order=2, length=2*pi) - out;

imshow(test);
colorbar();
figure();
imshow(out);
colorbar();
figure();
imshow(new-test);
colorbar();
show();


