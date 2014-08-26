from spectral import *;
from integrator import *;
from pylab import *;
import scipy as sp;
import pylab as m;


N = 256;
t = linspace(0,2*pi,N+1);
t = t[0:len(t)-1];
h = 2*pi/N;

tn = 0;
sec = 1;

#v = sin(t)+2*sin(t/4)+3*cos(t);
v = sin(t/2);
omega = zeros((N,N));
for i in arange(20,120):
    for j in arange(20,120):
        amplitude = 800*rand();
        phase = exp(j*rand());
        value = amplitude;
#        omega[N-i,N-j] = value;
        omega[i,j] = value;
        omega[N-i,j] = value;
        #omega[i,N-j] = value;

omega = real(ifft2(omega));

'''
omega = 3*rand(N,N);
omega = 3*zeros((N,N));
j = 0;
for i in arange(N):
    for j in arange(N):
        x = i-N/2;
        y = j-N/2;
        omega[i,j] = 5*exp(-(x*x + y*y)/10);
'''
'''
for i in arange(N):
    for j in arange(N/2):
        omega[i,j] += sin(4*pi*i/N);
'''
'''
for i in arange(N):
    for j in arange(N/2):
        omega[i,j+N/2] = rand();

'''
#omega += 4*rand(N,N);
'''
ensuring edges smoothly go to zero
for i in arange(N):
    omega[i,:] = omega[i,:]*v;

for i in arange(N):
    omega[:,i] = omega[:,i]*v;
'''

def dudt(dt, omega, args):
    """ 
        (psi_xx + psi_yy) = omega. Invert, calculate u = -psi_y, v = psi_x
        omega_t = -(u*omega_x + v*omega_y)
    """
    psi = InvLaplacian(omega);
    u = -partialY(psi);
    v = partialX(psi);

    return real(-(u*partialX(omega) + v*partialY(omega)));

def diffusion(dt, omega, args):
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

 #   delsq = -(KX*KX + KY*KY)
    delsq = -(pow(KX,6) + pow(KY,6))
    omegasq = sqrt(sum(omega*omega)/(N*N));
    print sum(omega*omega);

    nu = omegasq/pow(N,4);
    print 'viscosity=',nu;
    
    delsq = nu*delsq*dt;
    multiplier = exp(delsq);
    omegaHat = omegaHat*multiplier;
    return real(ifft2(omegaHat));

def calcMax(omega, args):
    psi = InvLaplacian(omega);
    u = -partialY(psi);
    v = partialX(psi);
    maximum = sqrt(amax(amax(u*u),amax(v*v)));
    print 'Max vel: ', maximum;
    return maximum;

def calcPower(field):
    [x,y] = shape(field);
    power = zeros(y);

    for i in arange(x):
        power = power + abs(fft(field[i,:]));

    power = power/x;
    return log(power);

stepper = integrator(h, [dudt], [diffusion], calcMax);
ion();
k = linspace(1,128,128);
omega_n = omega;
figure(figsize= (8,8));
subplot(211);
handle, = plot(log(k),calcPower(omega_n)[0:N/2]);
hold(True);
ylim([-8,8]);

# making my own colormap
cdict = {
    'red'  :  ((0., 0., 0.), (0.5, 0.5, 0.5), (1., 1, 1)),
    'green':  ((0., 0., 0.), (0.5, 0.5, 0.5), (1., 1, 1)),
    'blue' :  ((0., 0., 0.), (0.5, 0.5, 0.5), (1., 1., 1.))
}
#cdict = {
#    'red'  :  ((0., 0., 0.), (0.5, 0.25, 0.25), (1., 1., 1.)),
#    'green':  ((0., 1., 1.), (0.7, 0.0, 0.5), (1., 1., 1.)),
#    'blue' :  ((0., 1., 1.), (0.5, 0.0, 0.0), (1., 1., 1.))
#}

my_cmap = m.matplotlib.colors.LinearSegmentedColormap('my_colormap', cdict, 1024)



subplot(212);
handle2 = imshow(omega,vmin=-0.1,vmax=0.1,cmap=my_cmap);
tick = linspace(-10,10,20);
c2 = colorbar();
draw();



#jet();
while (sec < 3000):
    (tn,omega_n) = stepper.integrate(tn,omega_n);
    print 'Time elapsed: ',tn;
    if(tn > sec):
        handle.set_ydata(calcPower(omega_n)[0:N/2]);
       # handle.autoscale();
       # c1.update_bruteforce(handle);
        handle2.set_array(omega_n);
        #handle2.autoscale();
        c2.update_bruteforce(handle2);
        sec = sec+1;
        savefig('/home/joymm/figures/fig'+str(sec)+'.png');
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
