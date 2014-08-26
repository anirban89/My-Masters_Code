from spectral import *;
from integrator import *;
from pylab import *;
import scipy as sp;
import pylab as m;
import sys;

resume = 0;
path = '/home/joy/EqVaportestSink';

if(len(sys.argv) > 1):
    if(len(sys.argv) == 3):
    
        print "Resuming from";
        resume = int(sys.argv[2]);
        print resume;
        coeff = float(sys.argv[1]);
    if(len(sys.argv) == 2):
    
        coeff = float(sys.argv[1]);
        print "Kappa = ", str(coeff);
else:
    print 'Insufficient arguments';
    exit();

N = 128;
h = 2*pi/N;

tn = 0;
sec = 1;

'''max eta that should be reached'''
eta_u = 1;

saturatedField = zeros((N,N));
#v = sin(t)+2*sin(t/4)+3*cos(t);
omega = ones((N,N));
vapor = zeros((N,N));
#for i in arange(20,80):
#    for j in arange(20,80):
#        omega[i,j] = 1000*rand();
#        omega[N-i,j] = value;
#        omega[i,N-j] = value;
x = linspace(0,2*pi,N+1);
x = x[0:len(x)-1];
y = linspace(0,2*pi,N+1);
y = y[0:len(y)-1];



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
    (u,v,delta) = args;

    return real(-(u*partialX(omega) + v*partialY(omega) + v - delta));
    #return real(-(u*partialX(omega) + v*partialY(omega) + v));
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
    downfilter = exp((arange(numvals/4+1, (numvals/2)+1) -
    numvals/2.0))[::-1];
    upfilter = exp(-(arange(numvals/2+1, (3*numvals/4)+1) -
    numvals/2.0))[::-1];
    #downfilter = slopeDown*(arange(numvals/6+1, (numvals/2)+1) - numvals/2.0);
    #upfilter = slopeUp*(arange(numvals/2, (5*numvals/6)+1) - numvals/2.0);

    preserveVals = ones(numvals/4);

    fftFilter = concatenate((preserveVals, downfilter, upfilter, preserveVals))

    omegaHat = (fftFilter * omegaHat.transpose()).transpose();
    omegaHat = (fftFilter * omegaHat);


#    omegaHat[xVal/3:2*xVal/3,yVal/3:2*yVal/3] = 0;

    k = array(range(N), dtype=complex128) - N/2;
    k = ifftshift(k);
    [KX, KY] = meshgrid(k,k);

    mu = 0.005;
    mask = zeros(N);
    for i in arange(int(N/100)):
        mask[i] = 1;
        mask[N-i-1] = 1;

    [dragx,dragy] = meshgrid(mask,mask);


#   delsq = -(KX*KX + KY*KY)
    delsq = -(pow(KX,6) + pow(KY,6))
    omegasq = sqrt(sum(omega*omega)/(N*N));
#    print sum(omega*omega);

    nu = 10*omegasq/pow(N,4);
    print 'viscosity=',nu;
    
    delsq = (nu*delsq - mu*(dragx+dragy))*dt;
    multiplier = exp(delsq);
    omegaHat = omegaHat*multiplier;
    return real(ifft2(omegaHat));

def calcMax_omega(fields, args):
    (omega,vapor) = fields;
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
    (u,v,kappa,delta) = args;

    return real(-(u*partialX(vapor) + v*partialY(vapor) - kappa*v + delta));
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
#    return vapor;
    vaporHat = fft2(vapor);
    [xVal,yVal] = shape(vaporHat);
    numvals = len(vapor);
    downfilter = exp((arange(numvals/4+1, (numvals/2)+1) -
    numvals/2.0))[::-1];
    upfilter = exp(-(arange(numvals/2+1, (3*numvals/4)+1) -
    numvals/2.0))[::-1];
    #downfilter = slopeDown*(arange(numvals/6+1, (numvals/2)+1) - numvals/2.0);
    #upfilter = slopeUp*(arange(numvals/2, (5*numvals/6)+1) - numvals/2.0);

    preserveVals = ones(numvals/4);

    fftFilter = concatenate((preserveVals, downfilter, upfilter, preserveVals))

    vaporHat = (fftFilter * vaporHat.transpose()).transpose();
    vaporHat = (fftFilter * vaporHat);


#    vaporHat[xVal/3:2*xVal/3,yVal/3:2*yVal/3] = 0;

    k = array(range(N), dtype=complex128) - N/2;
    k = ifftshift(k);
    [KX, KY] = meshgrid(k,k);

    mu = 0.005;
    mask = zeros(N);
    for i in arange(int(N/100)):
        mask[i] = 1;
        mask[N-i-1] = 1;

    [dragx,dragy] = meshgrid(mask,mask);


#   delsq = -(KX*KX + KY*KY)
    delsq = -(pow(KX,6) + pow(KY,6))
    vaporsq = sqrt(sum(vapor*vapor)/(N*N));
#    print sum(vapor*vapor);

    nu = 40*vaporsq/pow(N,4);
    print 'viscosity vapor=',nu;
    
    delsq = (nu*delsq - mu*(dragx+dragy))*dt;
    multiplier = exp(delsq);
    vaporHat = vaporHat*multiplier;
    return real(ifft2(vaporHat));

def dSystemdt(dt,fields,args):
    [omega_new, vapor_new] = fields;
    (coeff2,forcing) = args;
    psi = InvPotentialVorticity(omega_new);
    u = -partialY(psi);
    v = partialX(psi);

#    delta = 1/(1+exp(-1000*vapor_new));
    delta = zeros((N,N));
    delta[vapor_new>0] = 1;
    delta = delta*vapor_new;

    forcingField = delta-forcing;

    dOmega = dudt_omega(dt, omega_new, (u,v,delta));
    dVapor = dudt_vapor(dt, vapor_new, (u,v,coeff2,forcingField));
    return((dOmega, dVapor));


def diffusion(dt, fields, args):
    [omega_new, vapor_new] = fields;

    diffOmega = diffusion_omega(dt, omega_new, args);
    diffVapor = diffusion_vapor(dt, vapor_new, args);

    return([diffOmega, diffVapor]);


stepper = integrator(h, [dSystemdt], \
                [diffusion], calcMax_omega,2);

ion();
k = linspace(1,128,128);
vapor_n = vapor;
omega_n = omega;

dy = 2*pi/(N-1);
qs = zeros((N));
q0 = 7;

betaMax = q0/(2*pi);

beta = coeff*betaMax;

for i in arange(N):
    qs[i] = q0 - beta*dy*(i);



for i in arange(N):
    saturatedField[:,i] = qs[i];

print qs;
#for i in arange(N):

#   for j in arange(N):
#        vapor_n[i,j] = -qs[j]*(rand());
#         vapor_n[i,j] = -0.1*rand();
    #    vapor_n[i,j] = exp((-pow(i-N/2,2)-pow(j-N/2,2))/200);
        
        #vapor_n[i,j] = rand()-0.8;

#vapor_n = vapor_n - 0.5;

#for j in arange(N/2-10,N/2+10):
#    for i in arange(N):

#        vapor_n[i,j] = 0;

#for i in arange(N):
#    for j in arange(N):
#        if(vapor_n[i,j] < -(qs[j])):
#            vapor_n[i,j] = -qs[j];

#        rh[i,j] = (vapor_n[i,j]/qs[j]) + 1;
#        if (rh[i,j] > 1):
#            rh[i,j] = 1;


figure(figsize= (16,8));

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
tn2 = 0;
sec = 0;
#jet();

if(resume > 0):
   print 'Resuming from step: ',resume;
   omega_n = \
       np.load(path+str(coeff)+'/eta_data_'+str(resume)+'.npy').transpose();
   vapor_n = \
       np.load(path+str(coeff)+'/vapor_data_'+str(resume)+'.npy').transpose();
   sec = resume+10;
   tn = resume;


subplot(121);
handle2 = imshow(omega_n.transpose(),vmin=-0.1,vmax=0.1);
tick = linspace(-10,10,20);
c2 = colorbar();
draw();

subplot(122);
handle3 = imshow(vapor_n.transpose(), vmin=-15,vmax=1);
#handle3 = imshow(vapor_n.transpose(), vmin=-15,vmax=1,cmap=my_cmap_vapor);
tick = linspace(-10,10,20);
c3 = colorbar();
draw();

forcingTime = 0;
global_time = 0;

while (sec<100000):
    psi = InvPotentialVorticity(omega_n);
    '''ensure that psi (eta) does not cross the ranges under which shallow water
    equations are valid'''

    #psi[psi>eta_u] = eta_u;
    #psi[psi<-eta_u] = -eta_u;

    u = -partialY(psi);
    v = partialX(psi);
    condensationCount = 0;
    negativeRHCount = 0;
    forcing = rand(N,N);
    (dt,fields) = stepper.integrate(tn,[omega_n, vapor_n],(coeff,forcing));

    [omega_n,vapor_n] = fields;

    print 'Max vapor: ', amax(amax(vapor_n));
    for i in arange(N):
        for j in arange(N):
#            if(vapor_n[i,j] > 0):
#                condensationCount = condensationCount+1;
                #delta[i,j] = 10*vapor_n[i,j]/qs[j];
#                delta[i,j] = 1;
                #vapor_n[i,j] = 0;
#            else:
#                delta[i,j] = 0;

            if(vapor_n[i,j] < -(qs[j])):
                vapor_n[i,j] = -qs[j];
                negativeRHCount = negativeRHCount+1;
                ''' hack to ensure that RH is non-negative always.
                ''' 

    tn = dt + prev_t;

#    delta[etaDiff > eta_u] = 0;

    ''' condensation strength is proportional to eta height'''

    #delta_eta = delta*etaDiff;
    '''
    if(tn < 3):
        delta_eta = delta*vapor_n;
    else:
        delta_eta = zeros((N,N));
    '''
#    print 'minimum vapor value: ', amin(amin(vapor));

    if(tn > sec):
        forcingTime = forcingTime+1;

        print '==================================================================='

        print 'Time elapsed: ',tn;
        print 'negative RH cells: ', negativeRHCount;


        print '==================================================================='


 #       for i in arange(N):
 #               for j in arange(N):
 #                   val = rand()-0.995;
 #                   if(val > 0):
 #                       vapor_n[i][j] = 0;
       # handle.autoscale();
       # c1.update_bruteforce(handle);
        height = InvPotentialVorticity(omega_n);
        np.save(path+str(coeff)+'/vapor_data_'+str(sec)+'.npy',vapor_n.transpose());
        np.save(path+str(coeff)+'/eta_data_'+str(sec)+'.npy',omega_n.transpose());
        #savefig(path+str(coeff)+'/fig'+str(sec)+'.png');
        handle3.set_array(vapor_n.transpose());
        handle3.autoscale();
        handle2.set_array(height.transpose());
        handle2.autoscale();
        c2.update_bruteforce(handle2);
        sec = sec+10;
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
