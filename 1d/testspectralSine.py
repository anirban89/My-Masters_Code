from pylab import *;
from spectral import *;

N = 2048;
length = 8*pi;
t = linspace(0,length,N+2);
t = t[1:len(t)-1];

print 'sin', shape(sin(t)), 'test', shape(test);
test = zeros((1,N));
for i in arange(N):
    test[:,i] = sin(t[i]);


out1 = partialYSine(test,1,length).transpose();
out2 = partialYSine(test,2,length).transpose();
out3 = partialYSine(test,3,length).transpose();
out4 = partialYSine(test,4,length).transpose();

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
