from pylab import *;
from spectral import *;

N = 2048;
length = 8*pi;
t = linspace(0,length,N+1);
t = t[0:len(t)-1];

print 'sin', shape(sin(t)), 'test', shape(test);
test = cos(t);


out1 = partialX(test,1,length);
out2 = partialX(test,2,length);
out3 = partialX(test,3,length);
out4 = partialX(test,4,length);

print amax(test)/amax(out1);
print amax(test)/amax(out2);
print amax(test)/amax(out3);
print amax(test)/amax(out4);
hold(True);
subplot(521);
plot(test);
subplot(522);
plot(abs(fourier(test)));
subplot(523);
plot(out1);
subplot(524);
plot(abs(fourier(out1)));
subplot(525);
plot(out2);
subplot(526);
plot(abs(fourier(out2)));
subplot(527);
plot(out3);
subplot(528);
plot(abs(fourier(out3)));
subplot(529);
plot(out4);
show();
