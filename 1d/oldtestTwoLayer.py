from spectral import *;
from integrator import *;
from pylab import *;
import scipy as sp;
import pylab as m;
import sys;

resume = 0;
path = '/home/joy/TwoLayer';

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

if (coeff < 0):
    path = '/home/joy/TwoLayerNegative';
N = 256;
h = 2*pi/N;

tn = 0;
sec = 1;

path = path+str(abs(coeff));

'''max eta that should be reached'''
eta_u = 1;

saturatedField = zeros((N,N));
#v = sin(t)+2*sin(t/4)+3*cos(t);
omega1 = 0.01*rand(N,N);
omega2 = 0.01*rand(N,N);
vapor = ones((N,N));
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
#    omega1[i,:] = 0.1*cos(7*x[i])*sin(7*y);
#    omega2[i,:] = 0.1*sin(7*x[i])*cos(7*y);


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

def domega1dt(dt, omega, args):
    """ 
        (psi_xx + psi_yy - psi) = omega. Invert, calculate u = -psi_y, v = psi_x
        omega_t = -(u*omega_x + v*omega_y + v + delta*psi)
    """
    (u,v,delta) = args;

    return real(-(u*partialX(omega) + v*partialY(omega) + v + delta));

def domega2dt(dt, omega, args):
    """ 
        (psi_xx + psi_yy - psi) = omega. Invert, calculate u = -psi_y, v = psi_x
        omega_t = -(u*omega_x + v*omega_y + v + delta*psi)
    """
    (u,v,delta) = args;

    return real(-(u*partialX(omega) + v*partialY(omega) + v - delta));

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

    mu = 0.005;
    mask = zeros(N);
    for i in arange(int(N/100)):
        mask[i] = 1;
        mask[N-i-1] = 1;

    [dragx,dragy] = meshgrid(mask,mask);


    delsq = -(pow(KX,6) + pow(KY,6))
    omegasq = sqrt(sum(omega*omega)/(N*N));

    nu = 10*omegasq/pow(N,4);
    print 'viscosity=',nu;
    
    delsq = (nu*delsq - mu*(dragx+dragy))*dt;
    multiplier = exp(delsq);
    omegaHat = omegaHat*multiplier;
    return real(ifft2(omegaHat));

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

def dudt_vapor(dt, vapor, args):
    """ 
        vapor_t = -(u*vapor_x + v*vapor_y + Beta*v)
    """
    (u,v,kappa,delta) = args;

    return real(-(u*partialX(vapor) + v*partialY(vapor) + kappa*u + delta));
    #return real(-( kappa*u + delta));
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

def deriveVelocities(omega1, omega2):
    '''
    omega1 = laplacian(psi1) + (psi2-psi1)
    omega2 = laplacian(psi2) + (psi1-psi2)
    '''
    sumPsi = InvLaplacian(omega1+omega2);
    diffPsi = (omega2 - omega1)/2.0;

    psi1 = (sumPsi + diffPsi)/2.0;
    psi2 = (sumPsi - diffPsi)/2.0;

    return(psi1,psi2);

def dSystemdt(dt,fields,args):
    [omega1_new, omega2_new, vapor_new] = fields;

    (coeff2,forcing) = args;
    (psi1, psi2) = deriveVelocities(omega1_new, omega2_new);
    u1 = -partialY(psi1);
    v1 = partialX(psi1);

    u2 = -partialY(psi2);
    v2 = partialX(psi2);

#    delta = 1/(1+exp(-1000*vapor_new));
    delta = zeros((N,N));
    delta[vapor_new>0] = 1;
    delta = delta*vapor_new;

    forcingField = delta-forcing;

    dOmega1 = domega1dt(dt, omega1_new, (u1,v1,delta));
    dOmega2 = domega2dt(dt, omega2_new, (u2,v2,delta));
    dVapor = dudt_vapor(dt, vapor_new, (u2,v2,coeff2,forcingField));
    return((dOmega1, dOmega2, dVapor));


def diffusion(dt, fields, args):
    [omega1_new, omega2_new, vapor_new] = fields;

    diffOmega1 = diffusion_omega(dt, omega1_new, args);
    diffOmega2 = diffusion_omega(dt, omega2_new, args);
    diffVapor = diffusion_vapor(dt, vapor_new, args);

    return([diffOmega1, diffOmega2, diffVapor]);


stepper = integrator(h, [dSystemdt], \
                [diffusion], calcMax_omega,3);

ion();
vapor_n = vapor;
omega1_n = omega1;
omega2_n = omega2;

dx = 2*pi/(N-1);
qs = zeros((N));
q0 = 7;

betaMax = q0/(2*pi);

beta = coeff*betaMax;

if(coeff > 0):
    for i in arange(N):
        qs[i] = q0 - beta*dx*(i);
elif(coeff < 0):
    for i in arange(N):
        qs[i] = -beta*dx*i;


for i in arange(N):
    saturatedField[i,:] = qs[i];

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


figure(figsize= (10,8));

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


subplot(111);
handle2 = imshow(vapor_n.transpose(), vmin=-15,vmax=1);
#handle3 = imshow(vapor_n.transpose(), vmin=-15,vmax=1,cmap=my_cmap_vapor);
tick = linspace(-10,10,20);
c3 = colorbar();
draw();

forcingTime = 0;
global_time = 0;

while (sec<1000):

    condensationCount = 0;
    negativeRHCount = 0;
    #forcing = 0.03*rand(N,N);
    forcing = zeros((N,N));
    (dt,fields) = stepper.integrate(tn,[omega1_n, omega2_n, vapor_n],(coeff,forcing));

    [omega1_n, omega2_n,vapor_n] = fields;

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

            if(vapor_n[j,i] < -(qs[j])):
                vapor_n[j,i] = -qs[j];
                negativeRHCount = negativeRHCount+1;
                ''' hack to ensure that RH is non-negative always.
                ''' 

    tn = dt + prev_t;
    prev_t = tn;
    print '============================='
    print 'Time Elapsed: ', tn;
    print '============================='

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
        np.save(path+'/vapor_data_'+str(sec)+'.npy',vapor_n.transpose());
        np.save(path+'/eta1_data_'+str(sec)+'.npy',omega1_n.transpose());
        np.save(path+'/eta2_data_'+str(sec)+'.npy',omega2_n.transpose());
        savefig(path+'/fig'+str(sec)+'.png');
        #handle1.set_array(omega1_n.transpose());
        #handle1.autoscale();
        #handle3.set_array(omega2_n.transpose());
        #handle3.autoscale();
        handle2.set_array(vapor_n.transpose());
        handle2.autoscale();
        #c2.update_bruteforce(handle2);
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
