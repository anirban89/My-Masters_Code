from spectral import *;
from pylab import *;

N = 1280;
t = linspace(0,2*pi,N+1);
t = t[0:len(t)-1];

test = zeros((N,N));
print 'sin', shape(sin(t)), 'test', shape(test);

for i in arange(N):
    test[i,:] = sin(15*t[i])*cos(15*t);


out = InvLaplacian(test,length=2*pi);

new = partialY(out,order=2,length=2*pi) + \
    partialX(out,order=2, length=2*pi);

imshow(test);
colorbar();
figure();
imshow(out);
colorbar();
figure();
imshow(new-test);
colorbar();
show();


