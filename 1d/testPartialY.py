from pylab import *
from finiteDifference import *
from boundary import *

N = 10;

a = ones((N,N));

north = zeros((N,4));
south = zeros((N,4));

bc = boundary(yBcType='Dirichlet',bcs=(north,south));

print boydPartialY(a,bc);
