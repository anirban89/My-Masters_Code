from spectral import *;
from diffusion import spectralDiffusion;
from integrator import *;
from pylab import *;
import scipy as sp;
import pylab as m;
import sys;
from os.path import expanduser;


home = expanduser("~");

resume = 0;
path = home+'/TwoLayerLinear';

if(len(sys.argv) > 1):
    if(len(sys.argv) == 4):
    
        print "Resuming from";
        resume = int(sys.argv[3]);
        print resume;
        alphaVapor = float(sys.argv[1]);
        betaVapor = float(sys.argv[2]);
        print "alphaQ = ", alphaVapor;
        print "betaQ = ", betaVapor;
    if(len(sys.argv) == 3):
    
        alphaVapor = float(sys.argv[1]);
        betaVapor = float(sys.argv[2]);
        print "alphaQ = ", alphaVapor;
        print "betaQ = ", betaVapor;
else:
    print 'Insufficient arguments';
    quit();

if (alphaVapor < 0):
    path = home+'/TwoLayerLinearNegative';
N = 256;
h = 2*pi/N;

tn = 0;
sec = 1;

path = path+str(int(10*abs(alphaVapor)))+str(int(10*betaVapor));

'''max eta that should be reached'''
eta_u = 1;

#v = sin(t)+2*sin(t/4)+3*cos(t);

Uimp = -0.0;

psi1 = 0.0*rand(N,N);
psi2 = 0.0*rand(N,N);
vapor = 3*rand(N,N);
slope = 0.1;

#for i in arange(N):
#    for j in arange(N):
#        vapor[i,j] = \
#                 3*np.random.rand(N,N);
                 #4+(sin(0*i*h)*sin(0*j*h)+sin(7*i*h)*sin(7*j*h)+sin(10*i*h)*sin(10*j*h));
                 #3*exp(-(pow((i-N/2)*h,2)/0.1+pow((j-N/2)*h,2)))*cos(3*i*h);
 #                0.1*exp(-pow((j-N/2)*h,2)/0.25);
 #              0.01*cos(2*pi*(N/2-i)*2/float(N))*sin(2*pi*(j)*0.5/float(N));
#                0.05*cos(2*pi*(N/2-i)*5/float(N));
psi1 = zeros((N,N));
vapor[vapor<0] *= -1;



omega1 = partialX(psi1,2) + partialY(psi1,2) + (psi2-psi1);
omega2 = partialX(psi2,2) + partialY(psi2,2) + (psi1-psi2);
temp = zeros((N,N));


def domega1dt(dt, omega, args):
    """ 
        (psi_xx + psi_yy - psi) = omega. Invert, calculate u = -psi_y, v = psi_x
        omega_t = -(u*omega_x + v*omega_y + v + delta*psi)
    """
    (u,v,delta) = args;

    return real(-(v + delta));

def domega2dt(dt, omega, args):
    """ 
        (psi_xx + psi_yy - psi) = omega. Invert, calculate u = -psi_y, v = psi_x
        omega_t = -(u*omega_x + v*omega_y + v + delta*psi)
    """
    (u,v,delta) = args;

    return real(-(v - delta));

def dudt_vapor(dt, vapor, args):
    """ 
        vapor_t = -(u*vapor_x + v*vapor_y + Beta*v)
    """
    (u,v,alphaQ,betaQ,delta) = args;

    return real(-(-betaQ*v + alphaQ*u + delta));

def deriveVelocities(omega1, omega2):
    '''
    omega1 = laplacian(psi1) - (psi1-psi2)
    omega2 = laplacian(psi2) + (psi1-psi2)
    '''
    sumPsi = InvLaplacian(omega1+omega2);
    diffPsi = InvPotentialVorticityTwoLayer(omega1 - omega2);

    psi1 = (sumPsi + diffPsi)/2.0;
    psi2 = (sumPsi - diffPsi)/2.0;

    return(psi1,psi2);

def calcMax_omega(fields, args):
    (omega1_new, omega2_new,vapor) = fields;

    (psi1, psi2) = deriveVelocities(omega1_new, omega2_new);
    u1 = -partialY(psi1);
    v1 = partialX(psi1);

    u2 = -partialY(psi2);
    v2 = partialX(psi2);

    maximum = \
        sqrt(amax((amax(u1*u1),amax(v1*v1),amax(u2*u2),amax(v2*v2))));
    print 'Max vel: ', maximum;
    return maximum;

def dSystemdt(dt,fields,args):
    [omega1_new, omega2_new, vapor_new] = fields;

    (alphaQ,betaQ,forcing) = args;
    (psi1, psi2) = deriveVelocities(omega1_new, omega2_new);
    u1 = -partialY(psi1);
    v1 = partialX(psi1);

    u2 = -partialY(psi2);
    v2 = partialX(psi2);

    u2 = u2 + Uimp;

#    delta = 1/(1+exp(-1000*vapor_new));
    delta = zeros((N,N),dtype='double');
    delta[vapor_new>0] = 1.0;
    delta = delta*vapor_new;

    forcingField = delta-forcing;

    dOmega1 = domega1dt(dt, omega1_new, (u1,v1,delta));
    dOmega2 = domega2dt(dt, omega2_new, (u2,v2,delta));
    dVapor = dudt_vapor(dt, vapor_new, (u2,v2,alphaQ,betaQ,forcingField));
    return((dOmega1, dOmega2, dVapor));


def diffusion(dt, fields, args):
    [omega1_new, omega2_new, vapor_new] = fields;

    diffOmega1 = spectralDiffusion(dt, omega1_new, args);
    diffOmega2 = spectralDiffusion(dt, omega2_new, args);
    diffVapor = spectralDiffusion(dt, vapor_new, args);

    return([diffOmega1, diffOmega2, diffVapor]);


stepper = integrator(h, [dSystemdt], \
                [diffusion], calcMax_omega,3);

ion();
vapor_n = vapor;
omega1_n = omega1;
omega2_n = omega2;

dx = 2*pi/(N-1);
dy = 2*pi/(N-1);
qs = zeros((N,N));
q0 = 7;

betaMax = q0/(2*pi);
alphaMax = q0/(2*pi);

beta = betaVapor*betaMax;
alpha = (alphaVapor*alphaMax);

if(alphaVapor != 0):
    for i in arange(N):
        for j in arange(N):
            qs[i,j] = q0 - alpha*dx*(i) - beta*dy*j;
else:
    for i in arange(N):
        for j in arange(N):
            qs[i,j] = q0 - beta*dy*j;



figure(figsize= (20,8));

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
   omega1_n = \
       np.load(path+'/eta1_data_'+str(resume)+'.npy').transpose();
   omega2_n = \
       np.load(path+'/eta2_data_'+str(resume)+'.npy').transpose();
   vapor_n = \
       np.load(path+'/vapor_data_'+str(resume)+'.npy').transpose();
   sec = resume+1;
   prev_t = resume;


#subplot(131);
#handle1 = imshow(omega1_n.transpose(),vmin=-0.1,vmax=0.1);
#tick = linspace(-10,10,20);

#c2 = colorbar();
#draw();

#subplot(121);
#handle3 = imshow(omega2_n.transpose(),vmin=-0.1,vmax=0.1);
#tick = linspace(-10,10,20);
#c2 = colorbar();
#draw();


subplot(121);
handle1 = imshow(vapor_n.transpose(), vmin=-15,vmax=1, aspect='auto');
#handle3 = imshow(vapor_n.transpose(), vmin=-15,vmax=1,cmap=my_cmap_vapor);
tick = linspace(-10,10,20);
c3 = colorbar();
hold(True);
contour(qs.transpose());
subplot(122);
SF1,SF2 = deriveVelocities(omega1_n,omega2_n);
handle2 = imshow((SF1-SF2).transpose(),aspect='auto');
c4 = colorbar();
draw();

forcingTime = 0;
global_time = 0;

while (sec<500):

    condensationCount = 0;
    negativeRHCount = 0;
    #forcing = .3*rand(N,N);
    forcing = 0.0*ones((N,N));
    (dt,fields) = stepper.integrate(tn,[omega1_n, omega2_n, \
                                            vapor_n],(alphaVapor,betaVapor,forcing));

    [omega1_n, omega2_n,vapor_n] = fields;

    print 'Max vapor: ', amax(amax(vapor_n));
    for i in arange(N):
        for j in arange(N):

            if(vapor_n[j,i] < -(qs[j,i])):
                vapor_n[j,i] = -qs[j,i];
                negativeRHCount = negativeRHCount+1;
                ''' hack to ensure that RH is non-negative always.
                ''' 

    tn = dt + prev_t;
    prev_t = tn;
    print '============================='
    print 'Time Elapsed: ', tn;
    print '============================='


    if(tn > sec):
        forcingTime = forcingTime+1;

        print '==================================================================='

        print 'Time elapsed: ',tn;
        print 'negative RH cells: ', negativeRHCount;


        print '==================================================================='
        delta = zeros((N,N),dtype='double');
        delta[vapor_n>0] = 1.0;
        delta = delta*vapor_n;

        SF1,SF2 = deriveVelocities(omega1_n,omega2_n);
        v = partialX(SF2);
        u = -partialY(SF2);
        handle1.set_array(vapor_n.transpose());
        handle1.autoscale();
        handle2.set_array(-(-betaVapor*v + alphaVapor*u + delta).transpose());
        handle2.autoscale();

        sec = sec+1;
        draw();
        np.save(path+'/vapor_data_'+str(sec)+'.npy',vapor_n.transpose());
        np.save(path+'/eta1_data_'+str(sec)+'.npy',omega1_n.transpose());
        np.save(path+'/eta2_data_'+str(sec)+'.npy',omega2_n.transpose());
        savefig(path+'/fig'+str(sec)+'.png');

        ioff();
        """
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

    preserveVals = ones(numvals/4);

    fftFilter = concatenate((preserveVals, downfilter, upfilter, preserveVals))

    omegaHat = (fftFilter * omegaHat.transpose()).transpose();
    omegaHat = (fftFilter * omegaHat);



    k = array(range(N), dtype=complex128) - N/2;
    k = ifftshift(k);
    [KX, KY] = meshgrid(k,k);

    mu = 0.0005;
    mask = zeros(N);
    for i in arange(int(N/100)):
        mask[i] = 1;
        mask[N-i-1] = 1;

    [dragx,dragy] = meshgrid(mask,mask);


    delsq = -(pow(KX,6) + pow(KY,6))
    omegasq = sqrt(sum(omega*omega)/(N*N));

    nu = omegasq/pow(N,4);
    print 'viscosity=',nu;
    
    delsq = (nu*delsq - mu*(dragx+dragy))*dt;
    multiplier = exp(delsq);
    omegaHat = omegaHat*multiplier;
    return real(ifft2(omegaHat));
"""
"""
def diffusion_vapor(dt, vapor, args):
    '''
        d(omega)/dt = nu*((omega)_xx + (omega)_yy)
        => d(omega)Hat/(omega)Hat = -nu*(k_x^(2n) + k_y^(2n))*dt
        => (omega)Hat = e^(-nu*|k|^2n)dt*(omega)Hat_0

        Dynamic viscosity calculated according to the formula
        nu = omega_rms/k_max^(2n-2), where n is the viscosity power, 3 for
        biharmonic function, used here.
    '''
    vaporHat = fft2(vapor);
    [xVal,yVal] = shape(vaporHat);
    numvals = len(vapor);
    downfilter = exp((arange(numvals/4+1, (numvals/2)+1) -
    numvals/2.0))[::-1];
    upfilter = exp(-(arange(numvals/2+1, (3*numvals/4)+1) -
    numvals/2.0))[::-1];

    preserveVals = ones(numvals/4);

    fftFilter = concatenate((preserveVals, downfilter, upfilter, preserveVals))

    vaporHat = (fftFilter * vaporHat.transpose()).transpose();
    vaporHat = (fftFilter * vaporHat);



    k = array(range(N), dtype=complex128) - N/2;
    k = ifftshift(k);
    [KX, KY] = meshgrid(k,k);

    mu = 0.0005;
    mask = zeros(N);
    for i in arange(int(N/100)):
        mask[i] = 1;
        mask[N-i-1] = 1;

    [dragx,dragy] = meshgrid(mask,mask);


#   delsq = -(KX*KX + KY*KY)
    delsq = -(pow(KX,6) + pow(KY,6))
    vaporsq = sqrt(sum(vapor*vapor)/(N*N));
#    print sum(vapor*vapor);

    nu = vaporsq/pow(N,4);
    print 'viscosity vapor=',nu;
    
    delsq = (nu*delsq - mu*(dragx+dragy))*dt;
    multiplier = exp(delsq);
    vaporHat = vaporHat*multiplier;
    return real(ifft2(vaporHat));
"""
