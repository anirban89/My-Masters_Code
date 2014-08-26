from diffs2d import *;
from pylab import *;
from spectral import *;

N = 1280;
t = linspace(0,2*pi,N+1);
t = t[0:len(t)-1];

test = zeros((N,N));
print 'sin', shape(sin(t)), 'test', shape(test);
for i in arange(N):
    test[i,:] = cos(t);

diff = diffs2d(N);

out = diff.psderiv(test,order=2,length=2*pi,axis=1);
out2 = partial(test,order=2,length=2*pi,axis=1);

imshow(test);
colorbar();
figure();
imshow(out);
colorbar();
figure();
imshow(out2);
colorbar();
show();
