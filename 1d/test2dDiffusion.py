from spectral import *;
from integrator import *;
from pylab import *;
from diffusion import *;
import scipy as sp;


N = 128;
t = linspace(0,2*pi,N+1);
t = t[0:len(t)-1];
h = 2*pi/N;
v = zeros((N,N));

#v = sin(t)+2*sin(t/4)+3*cos(t);

def dudt(dt, u, args):
    """ du/dt = -u*u_x 
    """
#    return 0;
    return real(-u*partialX(u));


def calcMax(u, args):
    return amax((0.2,sqrt(amax(u*u))));

for i in arange(N):
    v[i,:] = 0.1*sin(t[i])*cos(t);

stepper = integrator(h, [dudt], [spectralDiffusion], calcMax);
ion();
vn = v;
subplot(211);
handle = imshow(abs(fft2(v)));
subplot(212);
handle2 = imshow(v.transpose());
show();

tn = 0;
sec = 1;
while (sec < 500):
    (tn,vn) = stepper.integrate(tn,vn);
    if(tn > sec):
        print 'Time elapsed: ', tn;
        handle.set_array(abs(fft2(vn)));
        handle2.set_array(vn.transpose());
        sec = sec+1;
        draw();

ioff();
'''figure();
plot((tvals),log(vvals));
[slope,intercept] =  polyfit((tvals),log(vvals),1);
hold(True);
plot((tvals),(slope*(tvals) + intercept));
show();
print log(tvals);
print log(vvals);
vals =  polyfit((tvals),log(vvals),1);
print 'coeffs=',  vals;
'''
