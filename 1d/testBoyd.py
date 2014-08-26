from finiteDifference import *;
from boundary import *;

N = 256;

test = ones((N,256));
diff = zeros((N,256));
north = zeros((N,4));
south = zeros((N,4));

bc = boundary('Periodic','Dirichlet',(north,south));

for i in arange(N):
    for j in arange(256):

        test[i,j] = cos(2*pi*i/(N))*sin(2*pi*j/(N));
        diff[i,j] = -sin(2*pi*i/(N))*sin(2*pi*j/(N));

out = boydPartialX(test);
outY = boydPartialY(test,bc);

imshow((out-diff));
colorbar();
figure();
imshow(outY);
colorbar();
show();
