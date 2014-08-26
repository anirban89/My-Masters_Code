from spectral import *;
from integrator import *;
from pylab import *;
import scipy as sp;
import pylab as m;


N = 512;
t = linspace(0,2*pi,N+1);
t = t[0:len(t)-1];
h = 2*pi/N;

tn = 0;
sec = 1;

'''max eta that should be reached'''
eta_u = 0.5;

#v = sin(t)+2*sin(t/4)+3*cos(t);
v = sin(t/2);
omega = zeros((N,N));
vapor = zeros((N,N));
for i in arange(20,80):
    for j in arange(20,80):
        omega[i,j] = 1000*rand();
#        omega[N-i,j] = value;
#        omega[i,N-j] = value;

for i in arange(N):
    for j in arange(N):
        vapor[i,j] = rand();
        #omega[i,j] = 2*rand();
vapor = vapor-0.5;
#omega = omega-1;
omega = real(ifft2(omega));


'''normalize eta'''
eta = InvPotentialVorticity(omega);
maxEta = amax(eta);
eta = eta*0.1/maxEta;
omega = partialX(partialX(eta) + partialY(partialY(eta) - eta;


'''
for i in arange(N):
    for j in arange(5):
        omega[i,j] = 0;
        omega[i,N-j-1] = 0;

for i in arange(5):
    for j in arange(N):
        omega[i,j] = 0;
        omega[N-i-1,j] = 0;
'''
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

def dudt_omega(dt, omega, args):
    """ 
        (psi_xx + psi_yy - psi) = omega. Invert, calculate u = -psi_y, v = psi_x
        omega_t = -(u*omega_x + v*omega_y + v + delta*psi)
    """
    delta = args[0];
    psi = InvPotentialVorticity(omega);
    u = -partialY(psi);
    v = partialX(psi);

    return real(-(u*partialX(omega) + v*partialY(omega) + v + delta*psi));


def diffusion_omega(dt, omega, args):
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
    b = 0;
    mask = zeros(N);
    for i in arange(N/4):
        mask[i] = 1;
        mask[N-i-1] = 1;

    [dragx,dragy] = meshgrid(mask,mask);
    a = zeros((N,N));


#   delsq = -(KX*KX + KY*KY)
    delsq = -(pow(KX,6) + pow(KY,6))
    omegasq = sqrt(sum(omega*omega)/(N*N));
    print sum(omega*omega);

    nu = omegasq/pow(N,4);
    print 'viscosity=',nu;
    
    delsq = (nu*delsq + b*(dragx + dragy))*dt;
    multiplier = exp(delsq);
    omegaHat = omegaHat*multiplier;
    return real(ifft2(omegaHat));

def calcMax_omega(omega, args):
    psi = InvPotentialVorticity(omega);
    u = -partialY(psi);
    v = partialX(psi);
    maximum = sqrt(amax(amax(u*u),amax(v*v)));
    print 'Max vel: ', maximum;
    return maximum;

def dudt_vapor(dt, vapor, args):
    """ 
        vapor_t = -(u*vapor_x + v*vapor_y + Beta*v)
    """
    (beta,u,v) = args;

    return real(-(u*partialX(omega) + v*partialY(omega) + beta*v));

def diffusion_vapor(dt, vapor, args):
    '''
        d(omega)/dt = nu*((omega)_xx + (omega)_yy)
        => d(omega)Hat/(omega)Hat = -nu*(k_x^(2n) + k_y^(2n))*dt
        => (omega)Hat = e^(-nu*|k|^2n)dt*(omega)Hat_0

        Dynamic viscosity calculated according to the formula
        nu = omega_rms/k_max^(2n-2), where n is the viscosity power, 3 for
        biharmonic function, used here.
    '''
    return vapor;
    vaporHat = fft2(vapor);
    [xVal,yVal] = shape(vaporHat);
    vaporHat[xVal/3:2*xVal/3,yVal/3:2*yVal/3] = 0;

    k = array(range(N), dtype=complex128) - N/2;
    k = ifftshift(k);
    [KX, KY] = meshgrid(k,k);

#   delsq = -(KX*KX + KY*KY)
    delsq = -(pow(KX,6) + pow(KY,6))
    vaporsq = sqrt(sum(vapor*vapor)/(N*N));

    nu = vaporsq/pow(N,4);
    
    delsq = nu*delsq*dt;
    multiplier = exp(delsq);
    vaporHat = vaporHat*multiplier;
    return real(ifft2(vaporHat));

def calcMax_vapor(vapor, args):
    (beta,u,v) = args;
    maximum = sqrt(amax(amax(u*u),amax(v*v)));
    return maximum;


def calcPower(field):
    [x,y] = shape(field);
    power = zeros(y);

    for i in arange(x):
        power = power + abs(fft(field[i,:]));

    power = power/x;
    return log(power);

stepper_omega = integrator(h, [dudt_omega], [diffusion_omega], calcMax_omega);
stepper_vapor = integrator(h, [dudt_vapor], [diffusion_vapor], calcMax_vapor);

ion();
k = linspace(1,128,128);
vapor_n = vapor;
omega_n = omega;
rh = zeros((N,N));
beta = 10;

dy = 2*pi/N;
qs = zeros((N,1));
q0 = 100;

for i in arange(N):
    qs[N-i-1] = beta*dy*i + q0;

for i in arange(N):
    for j in arange(N):
        vapor_n[i,j] = 0;


for i in arange(N):
    for j in arange(N):
        if(vapor_n[i,j] < -(qs[j])):
            vapor_n[i,j] = -qs[j];

        rh[i,j] = (vapor_n[i,j]/qs[j]) + 1;
        if (rh[i,j] > 1):
            rh[i,j] = 1;

(n,bins) = histogram(rh,bins=10,normed=True);

figure(figsize= (16,4));
(n,bins) = histogram(rh,bins=10,normed=True);
subplot(131);
handle, = plot(bins[0:len(bins)-1],n);
hold(True);
ylim([0,5]);

# making my own colormap
cdict = {
    'red'  :  ((0., 0., 0.), (0.5, 0.5, 0.5), (1., 1, 1)),
    'green':  ((0., 0., 0.), (0.5, 0.5, 0.5), (1., 1, 1)),
    'blue' :  ((0., 0., 0.), (0.5, 0.5, 0.5), (1., 1., 1.))
}
my_cmap = m.matplotlib.colors.LinearSegmentedColormap('my_colormap', cdict, 1024)

# making my own colormap
cdict_vapor = {
    'red'  :  ((0., 0., 0.), (0.5, 0.5, 0.5), (1., 1, 1)),
    'green':  ((0., 0., 0.), (0.5, 0.5, 0.5), (1., 1, 1)),
    'blue' :  ((0., 0., 0.), (0.5, 0.5, 0.5), (1, 0, 1))
}
my_cmap_vapor = m.matplotlib.colors.LinearSegmentedColormap('my_colormap', cdict_vapor, 1024)



subplot(132);
handle2 = imshow(omega.transpose(),vmin=-0.1,vmax=0.1,cmap=my_cmap);
tick = linspace(-10,10,20);
c2 = colorbar();
draw();

subplot(133);
handle3 = imshow(vapor_n.transpose(), vmin=-0.1,vmax=0.1,cmap=my_cmap_vapor);
tick = linspace(-10,10,20);
c3 = colorbar();
draw();

prev_t = 0;
delta = zeros((N,N));
tn2 = 0;
#jet();

while (sec < 3000):
    psi = InvPotentialVorticity(omega_n);
    u = -partialY(psi);
    v = partialX(psi);
    condensationCount = 0;
    negativeRHCount = 0;
    vorticity = partialX(v) - partialY(u);
    (tn,vapor_n) = stepper_vapor.integrate(tn,vapor_n,(beta,u,v));

    for i in arange(N):
        for j in arange(N):
            vapor[i,j] = vapor_n[i,j];
            delta[i,j] = 0;

    for i in arange(N):
        for j in arange(N):
            if(vapor_n[i,j] > 0):
                condensationCount = condensationCount+1;
                delta[i,j] = 1;
                vapor_n[i,j] = 0;
            else:
                delta[i,j] = 0;

            if(vapor_n[i,j] < -(qs[j])):
                vapor_n[i,j] = -qs[j];
                negativeRHCount = negativeRHCount+1;
                ''' hack to ensure that RH is non-negative always.
                ''' 

    dt = tn - prev_t;
    (tn2,omega_n) = stepper_omega.integrate(tn,omega_n,(delta));
    print  'v,q coherence: ', sum(v*vapor);
    print 'Time elapsed: ',tn;
    print 'Number of cells condensing: ', condensationCount;
    print 'negative RH cells: ', negativeRHCount;
    print 'minimum vapor value: ', amin(amin(vapor));

    if(tn > sec):
        rh = zeros((N,N));
        for i in arange(N):
            for j in arange(N):
                rh[i,j] = (vapor_n[i,j]/qs[j]) + 1;

        (n,bins) = histogram(rh,bins=10,normed=True);
        handle.set_ydata(n);
       # handle.autoscale();
       # c1.update_bruteforce(handle);
        np.save('/home/joymm/Vapor/vapor_data_'+str(sec)+'.npy',vapor.transpose());
        handle3.set_array(vapor.transpose());
        #handle3.autoscale();
        np.save('/home/joymm/Vapor/vorticity_data_'+str(sec)+'.npy',vorticity.transpose());
        handle2.set_array(vorticity.transpose());
        handle2.autoscale();
        c2.update_bruteforce(handle2);
        sec = sec+0.25;
        savefig('/home/joymm/Vapor/fig'+str(sec)+'.png');
        draw();

ioff();
logfile.close();
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
