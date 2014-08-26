from spectral import *;
from integrator import *;
from pylab import *;
import scipy as sp;


N = 1280;
t = linspace(0,2*pi,N+1);
t = t[0:len(t)-1];
h = 2*pi/N;
#v = sin(t)+2*sin(t/4)+3*cos(t);
v = sin(t);
#v[N/2:-1] = 0;

def dudt(dt, u, args):
    """ du/dt = -u*u_x 
    """
    return real(-u*partialX(u,1));

def diffusion(dt, u, args):
    '''
        du/dt = u_xx 
        => duHat/uHat = -k^2*dt
        => uHat = e^(-k^2)dt*uHat_0
    '''
    uHat = fourier(u);
    numvals = len(u);
    threshold = int(numvals/3);
    uHat[threshold:2*numvals/3] = 0;
    k = array(range(N), dtype=complex128) - N/2;
    k = ifftshift(k);
    k = 0.01*(k**2);
    k = -k*dt;
    multiplier = exp(k);
    uHat = uHat*multiplier;
    return real(invFourier(uHat));

def calcMax(u, args):
    return sqrt(amax(u*u));

stepper = integrator(h, [dudt], [diffusion], calcMax);
ion();
vn = v;
subplot(211);
handle, = plot(abs(fourier(v)));
subplot(212);
handle2, = plot(t,v);

tn = 0;
sec = 1;
while (sec < 30):
    (tn,vn) = stepper.integrate(tn,vn);
    if(tn > sec):
        handle.set_ydata((fourier(vn).real));
        handle2.set_ydata(vn);
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
