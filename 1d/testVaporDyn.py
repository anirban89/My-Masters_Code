from spectral import *;
from integrator import *;
from pylab import *;
import scipy as sp;
import pylab as m;
import sys;

resume = 0;
path = '/home/joy/Vapor';

if(len(sys.argv) > 1):
     print "Resuming";
     resume = int(sys.argv[1]);

N = 256;
t = linspace(0,2*pi,N+1);
t = t[0:len(t)-1];
h = 2*pi/N;

tn = 0;
sec = 1;

'''max eta that should be reached'''
eta_u = 1;

#v = sin(t)+2*sin(t/4)+3*cos(t);
v = sin(t/2);
omega = ones((N,N));
vapor = zeros((N,N));
#for i in arange(20,80):
#    for j in arange(20,80):
#        omega[i,j] = 1000*rand();
#        omega[N-i,j] = value;
#        omega[i,N-j] = value;


'''normalize eta'''
eta = InvPotentialVorticity(omega);
maxEta = amax(eta);
if(maxEta > 0):
    eta = eta*0.1/maxEta;
else:
    eta[:,:] = 0.00;

#for i in arange(N):
#    for j in arange(N):
#        omega[i,j] = exp(-(pow(i-N/2,2)+pow(j-N/2,2))/200);


#omega = partialX(partialX(eta)) + partialY(partialY(eta)) - eta;


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
    delta = args;
    psi = InvPotentialVorticity(omega);
    u = -partialY(psi);
    v = partialX(psi);

    return real(-(u*partialX(omega) + v*partialY(omega) + v + delta));
    #return real(-(u*partialX(omega) + v*partialY(omega) + v + delta*psi));


def diffusion_omega(dt, omega, args):
    '''
        d(omega)/dt = nu*((omega)_xx + (omega)_yy) - mu*omega
        => d(omega)Hat/(omega)Hat = -nu*(k_x^2 + k_y^2)*dt - mu*dt
        => (omega)Hat = e^(-nu*|k|^2 - mu)dt*(omega)Hat_0

        Dynamic viscosity nu = 0.01 for now
    '''
    omegaHat = fft2(omega);
    [xVal,yVal] = shape(omegaHat);
    numvals = len(omega);
    downfilter = exp((arange(numvals/6+1, (numvals/2)+1) -
    numvals/2.0))[::-1];
    upfilter = exp(-(arange(numvals/2, (5*numvals/6)+1) -
    numvals/2.0))[::-1];
    #downfilter = slopeDown*(arange(numvals/6+1, (numvals/2)+1) - numvals/2.0);
    #upfilter = slopeUp*(arange(numvals/2, (5*numvals/6)+1) - numvals/2.0);

    preserveVals = ones(numvals/6);

    fftFilter = concatenate((preserveVals, downfilter, upfilter, preserveVals))

    omegaHat = (fftFilter * omegaHat.transpose()).transpose();
    omegaHat = (fftFilter * omegaHat);


#    omegaHat[xVal/3:2*xVal/3,yVal/3:2*yVal/3] = 0;

    k = array(range(N), dtype=complex128) - N/2;
    k = ifftshift(k);
    [KX, KY] = meshgrid(k,k);

    mu = 0.0001;
    mask = zeros(N);
    for i in arange(int(N/100)):
        mask[i] = 1;
        mask[N-i-1] = 1;

    [dragx,dragy] = meshgrid(mask,mask);


#   delsq = -(KX*KX + KY*KY)
    delsq = -(pow(KX,6) + pow(KY,6))
    omegasq = sqrt(sum(omega*omega)/(N*N));
#    print sum(omega*omega);

    nu = omegasq/pow(N,4);
    print 'viscosity=',nu;
    
    delsq = (nu*delsq - mu*(dragx+dragy))*dt;
    multiplier = exp(delsq);
    omegaHat = omegaHat*multiplier;
    return real(ifft2(omegaHat));

def calcMax_omega(omega, args):
    psi = InvPotentialVorticity(omega);
    u = -partialY(psi);
    v = partialX(psi);
    maximum = sqrt(amax((amax(u*u),amax(v*v))));
    print 'Max vel: ', maximum;
    return maximum;

def dudt_vapor(dt, vapor, args):
    """ 
        vapor_t = -(u*vapor_x + v*vapor_y + Beta*v)
    """
    (beta,u,v,delta) = args;

    return real(-(u*partialX(vapor) + v*partialY(vapor) - beta*v + delta*vapor));
    #return real(-(u*partialX(vapor) + v*partialY(vapor) - beta*v));

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
    (beta,u,v,delta) = args;
    maximum = sqrt(amax((amax(u*u),amax(v*v))));
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

dy = 2*pi/N;
qs = zeros((N,1));
q0 = 7;

betaMax = q0/(2*pi);
coeff = 0.6;

beta = coeff*betaMax;

for i in arange(N):
    qs[i] = q0 - beta*dy*i;

for i in arange(N):
    for j in arange(N):
        vapor_n[i,j] = -qs[j]*(rand()-0.8);
        #vapor_n[i,j] = rand()-0.8;

#vapor_n = vapor_n - 0.5;

for i in arange(N):
    for j in arange(N):
        if(vapor_n[i,j] < -(qs[j])):
            vapor_n[i,j] = -qs[j];

        rh[i,j] = (vapor_n[i,j]/qs[j]) + 1;
        if (rh[i,j] > 1):
            rh[i,j] = 1;


figure(figsize= (16,4));
(n,bins) = \
histogram(rh,bins=[-0.05,0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95,1.05],normed=True);
subplot(131);
x = array([0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]);
handle, = plot(x,n,'o');
hold(True);
ylim([0,20]);

# making my own colormap
cdict = {
    'red'  :  ((0., 1., 1.), (0.5, 0.5, 0.5), (1., 0, 0)),
    'green':  ((0., 1., 1.), (0.5, 0.5, 0.5), (1., 0, 0)),
    'blue' :  ((0., 1., 1.), (0.5, 0.5, 0.5), (1., 0., 0.))
}
my_cmap = m.matplotlib.colors.LinearSegmentedColormap('my_colormap', cdict, 1024)

# making my own colormap
cdict_vapor = {
    'red'  :  ((0., 0, 0), (0.5, 0.5, 0.5), (1., 0.5, 1)),
    'green':  ((0., 1., 1.), (0.5, 0.5, 0.5), (1., 0, 0)),
    'blue' :  ((0., 0.5, 0), (0.5, 0.5, 0.5), (1, 0.5, 1))
}
my_cmap_vapor = m.matplotlib.colors.LinearSegmentedColormap('my_colormap', cdict_vapor, 1024)

prev_t = 0;
delta = zeros((N,N));
tn2 = 0;
sec = 0;
#jet();

if(resume > 0):
   print 'Resuming from step: ',resume;
   omega_n = \
       np.load(path+str(coeff)+'/eta_data_'+str(resume)+'.npy').transpose();
   vapor_n = \
       np.load(path+str(coeff)+'/vapor_data_'+str(resume)+'.npy').transpose();
   sec = resume+1;
   tn = resume;


subplot(132);
handle2 = imshow(omega_n.transpose(),vmin=-0.1,vmax=0.1,cmap=my_cmap);
tick = linspace(-10,10,20);
c2 = colorbar();
draw();

subplot(133);
handle3 = imshow(vapor_n.transpose(), vmin=-15,vmax=1);
#handle3 = imshow(vapor_n.transpose(), vmin=-15,vmax=1,cmap=my_cmap_vapor);
tick = linspace(-10,10,20);
c3 = colorbar();
draw();



while (sec<1500):
    psi = InvPotentialVorticity(omega_n);
    '''ensure that psi (eta) does not cross the ranges under which shallow water
    equations are valid'''

    #psi[psi>eta_u] = eta_u;
    #psi[psi<-eta_u] = -eta_u;

    u = -partialY(psi);
    v = partialX(psi);
    condensationCount = 0;
    negativeRHCount = 0;
    vorticity = partialX(v) - partialY(u);
    (tn,vapor_n) = stepper_vapor.integrate(tn,vapor_n,(coeff,u,v,delta));

    '''actual value to multiply in momentum equation if condensation occurs'''
    etaDiff = eta_u - psi;

    ''' if psi is greater than eta_u, then do not force system '''
    etaDiff[etaDiff<0] = 0;

    for i in arange(N):
        for j in arange(N):
            vapor[i,j] = vapor_n[i,j];
            delta[i,j] = 0;

    print 'Max vapor: ', amax(amax(vapor_n));
    for i in arange(N):
        for j in arange(N):
            if(vapor_n[i,j] > 0):
                condensationCount = condensationCount+1;
                #delta[i,j] = 10*vapor_n[i,j]/qs[j];
                delta[i,j] = 0.5;
                #vapor_n[i,j] = 0;
            else:
                delta[i,j] = 0;

            if(vapor_n[i,j] < -(qs[j])):
                vapor_n[i,j] = -qs[j];
                negativeRHCount = negativeRHCount+1;
                ''' hack to ensure that RH is non-negative always.
                ''' 

    dt = tn - prev_t;

#    delta[etaDiff > eta_u] = 0;

    ''' condensation strength is proportional to eta height'''

    #delta_eta = delta*etaDiff;
    delta_eta = delta*vapor_n;

    print 'forcing strength: ',sum(delta);

    (tn2,omega_n) = stepper_omega.integrate(tn,omega_n,(delta_eta));
    print  'v,q coherence: ', sum(v*vapor);
    print 'Time elapsed: ',tn;
    print 'Number of cells condensing: ', condensationCount;
    print 'negative RH cells: ', negativeRHCount;
#    print 'minimum vapor value: ', amin(amin(vapor));

    if(tn > sec):
        rh = zeros((N,N));
        for i in arange(N):
            for j in arange(N):
                rh[i,j] = (vapor_n[i,j]/qs[j]) + 1;

        (n,bins) = \
                   histogram(rh,bins=[-0.05,0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95,1.05],normed=True);
        handle.set_ydata(n);
       # handle.autoscale();
       # c1.update_bruteforce(handle);
        np.save('/home/joy/Vapor'+str(coeff)+'/rh_data_'+str(sec)+'.npy',rh.transpose());
        np.save('/home/joy/Vapor'+str(coeff)+'/vapor_data_'+str(sec)+'.npy',vapor_n.transpose());
        np.save('/home/joy/Vapor'+str(coeff)+'/delta_data_'+str(sec)+'.npy',delta.transpose());
        handle3.set_array(delta_eta.transpose());
        handle3.autoscale();
        np.save('/home/joy/Vapor'+str(coeff)+'/eta_data_'+str(sec)+'.npy',omega_n.transpose());
        handle2.set_array(omega_n.transpose());
        handle2.autoscale();
        c2.update_bruteforce(handle2);
        sec = sec+1;
        #savefig('/home/joy/Vapor'+str(coeff)+'/fig'+str(sec)+'.png');
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
