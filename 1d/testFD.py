from boundary import *;
from finiteDifference import *;
from pylab import *;

N = 500;
M = 256;

u = zeros((N-1,M-2));
uDiff = zeros((N-1,M-2));
nDiff = zeros((N-2,M-2));

for i in arange(N-1):
    for j in arange(M-2):

        u[i,j] = sin(2*pi*(i)/(N-1))*cos(2*pi*(j+0.5)/(M-2));
        uDiff[i,j] = -cos(2*pi*(i)/(N-1))*cos(2*pi*(j+0.5)/(M-2));

print 'Done initing u, uDiff';

for i in arange(N-2):
    for j in arange(M-2):
        nDiff[i,j] = -cos(2*pi*(i+0.5)/(N-1))*cos(2*pi*(j+0.5)/(M-2));

print 'Done initing nDiff';

bc = boundary();

uD = partialXfd(u, bc, 'u', 'n');

print shape(uD);
print shape(uDiff);

diff = nDiff - uD;

print diff

imshow(u.transpose());
colorbar();
figure();
imshow(uD.transpose());
colorbar();
figure();
imshow(nDiff.transpose());
colorbar();
figure();
imshow(diff.transpose());
colorbar();


show();
