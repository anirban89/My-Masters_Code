from spectral import *

N = 40;

t = linspace(0,pi,N);
x = cos(t)

kx,ky = meshgrid(x,x)

signal = sin(kx)
actual = cos(kx)

out = partialChebyshev(signal);


hold(True)

imshow(signal)
figure();
imshow(actual)
figure()
imshow(out)

print norm(abs(out-actual),inf);

show();
