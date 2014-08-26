from diffs1d import *;
from pylab import *;

N = 1280;
t = linspace(0,2*pi,N+1);
t = t[0:len(t)-1];

test = zeros((1,N));
print 'sin', shape(sin(t)), 'test', shape(test);
test = cos(t);

diff = diffs1d(N);

out = diff.partial(test);

hold(True);
plot(test);
plot(out);
show();
