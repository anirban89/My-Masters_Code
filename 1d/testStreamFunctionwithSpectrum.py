from spectral import *;
from integrator import *;
from pylab import *;
import scipy as sp;


N = 256;
t = linspace(0,2*pi,N+1);
t = t[0:len(t)-1];
h = 2*pi/N;
#v = sin(t)+2*sin(t/4)+3*cos(t);
v = sin(t/2);
omega = 3*rand(N,N);
omega = 3*zeros((N,N));

j = 0;
for i in arange(N):
    for j in arange(N):
        x = i-N/2;
        y = j-N/2;
        omega[i,j] = 2*exp(-(x*x + y*y)/10);

omega += rand(N,N);
'''
ensuring edges smoothly go to zero
for i in arange(N):
    omega[i,:] = omega[i,:]*v;

for i in arange(N):
    omega[:,i] = omega[:,i]*v;
'''

def dudt(dt, omega):
    """ 
        (psi_xx + psi_yy) = omega. Invert, calculate u = -psi_y, v = psi_x
        omega_t = -(u*omega_x + v*omega_y)
    """
    psi = InvLaplacian(omega);
    u = -partialY(psi);
    v = partialX(psi);

    return real(-(u*partialX(omega) + v*partialY(omega)));

def diffusion(dt, omega):
    '''
        d(omega)/dt = nu*((omega)_xx + (omega)_yy)
        => d(omega)Hat/(omega)Hat = -nu*(k_x^2 + k_y^2)*dt
        => (omega)Hat = e^(-nu*|k|^2)dt*(omega)Hat_0

        Dynamic viscosity nu = 0.01 for now
    '''
    omegaHat = fft2(omega);
    [xVal,yVal] = shape(omegaHat);
    omegaHat[xVal/3:2*xVal/3,yVal/3:2*yVal/3] = 0;

    k = array(range(N), dtype=complex128) - N/2;
    k = ifftshift(k);
    [KX, KY] = meshgrid(k,k);

    delsq = -(KX*KX + KY*KY)
    
    delsq = 0.00001*delsq*dt;
    multiplier = exp(delsq);
    omegaHat = omegaHat*multiplier;
    return real(ifft2(omegaHat));

def calcMax(omega):
    psi = InvLaplacian(omega);
    u = -partialY(psi);
    v = partialX(psi);
    return sqrt(amax(amax(u*u),amax(v*v)));

def calcPower(field):
    [x,y] = shape(field);
    power = zeros(y);

    for i in arange(x):
        power = power + abs(fft(field[i,:]));

    power = power/x;
    return power;

stepper = integrator(h, [dudt], [diffusion], calcMax);
ion();
omega_n = omega;
subplot(211);
handle, = plot(calcPower(omega_n));
subplot(212);
handle2 = imshow(omega);
tick = linspace(-10,10,20);
c2 = colorbar();
clim(0,0.5);
draw();
tn = 0;
sec = 1;
while (sec < 200):
    (tn,omega_n) = stepper.integrate(tn,omega_n);
    print 'Time elapsed: ',tn;
    if(tn > sec):
        handle.set_ydata(calcPower(omega_n));
       # handle.autoscale();
       # c1.update_bruteforce(handle);
        handle2.set_array(omega_n);
        handle2.autoscale();
        c2.update_bruteforce(handle2);
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
