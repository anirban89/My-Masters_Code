from spectral import *;
from integrator import *;
from pylab import *;
import scipy as sp;
import pylab as m;


N = 256;
t = linspace(0,pi,N/2+1);
t = t[1:len(t)-1];
h = 2*pi/N;

tn = 0;
sec = 1;

'''max eta that should be reached'''
eta_u = 1;

globalVar = 0;
#v = sin(t)+2*sin(t/4)+3*cos(t);
v = sin(t);
omega = ones((N,N/2+1));


k = linspace(1,128,128);
rh = zeros((N,N/2+1));
beta = 0.2;

dy = 2*pi/(N/2+1);
qs = zeros((N/2+1,1));
q0 = 0.3;

for i in arange(N/2+1):
    qs[N/2-i] = beta*dy*i + q0;



#for i in arange(N/2-1):
#    omega[:,i] = omega[:,i]*v[i];
vapor = zeros((N,N/2+1));
for i in arange(N):
    for j in arange(N/2+1):
        vapor[i,j] = qs[j]*(rand()+0.1);

#for i in arange(20,80):
#    for j in arange(20,80):
#        omega[i,j] = 1000*rand();
#        omega[N-i,j] = value;
#        omega[i,N-j] = value;


'''normalize eta'''
'''print 'NORMALIZING ETA...'
eta = InvPotentialVorticitySine(omega);
maxEta = amax(eta);
if(maxEta > 0):
    eta = eta*0.1/maxEta;
else:
    eta[:,:] = 0.00;
'''


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
    psi = InvPotentialVorticityCosine(omega);
    u = -partialYCosine(psi);
    v = partialX(psi);

    return real(-(u*partialX(omega) + v*partialYCosine(omega) + v + delta*psi));


def diffusion_omega(dt, omega, args):
    '''
        d(omega)/dt = nu*((omega)_xx + (omega)_yy) - mu*omega
        => d(omega)Hat/(omega)Hat = -nu*(k_x^2 + k_y^2)*dt - mu*dt
        => (omega)Hat = e^(-nu*|k|^2 - mu)dt*(omega)Hat_0

        Dynamic viscosity nu = 0.01 for now
    '''

    N = shape(omega)[0];
    leny = shape(omega)[1];
    newLen = 2*(leny-1);

    new = zeros((N,N));

    for i in arange(leny):
        new[:,i] = omega[:,i];

    for i in arange(leny-1):
        new[:,newLen-i-1] = omega[:,i+1];


    omegaHat = fft2(new);
    [xVal,yVal] = shape(omegaHat);
    omegaHat[xVal/3:2*xVal/3,yVal/3:2*yVal/3] = 0;

    k = array(range(N), dtype=complex128) - N/2;
    k = ifftshift(k);
    [KX, KY] = meshgrid(k,k);
    mu = 0.0001;
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
    
    delsq = (nu*delsq)*dt;
    multiplier = exp(delsq);
    omegaHat = omegaHat*multiplier;
    return real(ifft2(omegaHat)[:,0:leny]);

def calcMax_omega(omega, args):
    psi = InvPotentialVorticityCosine(omega);
    u = -partialYCosine(psi);
    v = partialX(psi);
    maximum = sqrt(amax(amax(u*u),amax(v*v)));
    print 'Max vel: ', maximum;
    return maximum;

var = 0;
def dudt_vapor(dt, vapor, args):
    """ 
        vapor_t = -(u*vapor_x + v*vapor_y + Beta*v)
    """
    (beta,u,v) = args;

    """ enforcing boundary conditions """
    for i in arange(15):
        vapor[:,i] = qs[i];

    global var;
#    np.save('/home/joymm/VaporS4/inter_data_'+str(var)+'.npy',vapor);
    var = var+1;

    return real(-(u*partialX(vapor) + v*partialYCosine(vapor)));

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
    N = shape(vapor)[0];
    leny = shape(vapor)[1];
    newLen = 2*(leny-1);

    for i in arange(10):
        vapor[:,i] = qs[i];
    new = zeros((N,N));

    for i in arange(leny):
        new[:,i] = vapor[:,i];
    for i in arange(leny-1):
        new[:,newLen-i-1] = vapor[:,i+1];


    vaporHat = fft2(new);
    [xVal,yVal] = shape(vaporHat);
    vaporHat[xVal/3:2*xVal/3,yVal/3:2*yVal/3] = 0;

    k = array(range(N), dtype=complex128) - N/2;
    k = ifftshift(k);
    [KX, KY] = meshgrid(k,k);

#   delsq = -(KX*KX + KY*KY)
    delsq = -(pow(KX,6) + pow(KY,6))
    vaporsq = sqrt(sum(vapor*vapor)/(N*N));

    nu = vaporsq/pow(N,4);
    
    delsq = 4*nu*delsq*dt;
    multiplier = exp(delsq);
    vaporHat = vaporHat*multiplier;
    return real(ifft2(vaporHat)[:,1:leny+1]);

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
vapor_n = vapor;
omega_n = omega;
#for j in arange(N/2+1):
#    vapor_n[:,j] = qs[j];

#for i in arange(N/4,3*N/4):
#    for j in arange(N/8,3*N/8):
#        vapor_n[i,j] = exp(-(pow(i-N/2,2)+pow(j-N/4-2,2))/200);

#for i in arange(N):
#    for j in arange(N/2-1):
#        vapor_n[i,j] = vapor[i,j]*qs[0];

for i in arange(N):
    for j in arange(N/2+1):
        if(vapor_n[i,j] < 0):
            vapor_n[i,j] = 0;

        rh[i,j] = (vapor_n[i,j]/qs[j]);
        if (rh[i,j] > 1):
            rh[i,j] = 1;


figure(figsize= (16,8));
(n,bins) = \
histogram(rh,bins=[-0.05,0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95,1.05],normed=True);
subplot(231);
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



subplot(232);
handle2 = imshow(omega.transpose(),vmin=-0.1,vmax=0.1,cmap=my_cmap);
tick = linspace(-10,10,20);
c2 = colorbar();
draw();

subplot(233);
handle3 = imshow(vapor_n.transpose(), vmin=-15,vmax=1,cmap=my_cmap_vapor);
tick = linspace(-10,10,20);
c3 = colorbar();
draw();

sfn = InvPotentialVorticityCosine(omega);
u = -partialYCosine(omega);
v = partialX(omega);

subplot(234);
handle4 = imshow(u.transpose(), vmin=-15,vmax=1,cmap=my_cmap_vapor);
tick = linspace(-10,10,20);
c3 = colorbar();
draw();

subplot(235);
handle5 = imshow(u.transpose(), vmin=-15,vmax=1,cmap=my_cmap_vapor);
tick = linspace(-10,10,20);
c3 = colorbar();
draw();

prev_t = 0;
delta = zeros((N,N/2+1));
tn2 = 0;
#jet();



np.save('/home/joymm/VaporS4/eta_data_0.npy',omega_n.transpose());
#while (sec < 1000):
while (1):
    psi = InvPotentialVorticityCosine(omega_n);
    '''ensure that psi (eta) does not cross the ranges under which shallow water
    equations are valid'''

   # psi[psi>eta_u] = eta_u;
   # psi[psi<-eta_u] = -eta_u;

    ''' Force boundary conditions: bottom always saturated'''

    for i in arange(15):
        vapor_n[:,i] = qs[i];



    u = -partialYCosine(psi);
    v = partialX(psi);
    condensationCount = 0;
    negativeRHCount = 0;
    vorticity = partialX(v) - partialYCosine(u);
    (tn,vapor_n) = stepper_vapor.integrate(tn,vapor_n,(beta,u,v));


    '''actual value to multiply in momentum equation if condensation occurs'''
    etaDiff = eta_u - psi;

    ''' if psi is greater than eta_u, then do not force system '''
    etaDiff[etaDiff<0] = 0;

    for i in arange(N):
        for j in arange(N/2+1):
            vapor[i,j] = vapor_n[i,j];
            delta[i,j] = 0;

    for i in arange(N):
        for j in arange(N/2+1):
            if(vapor_n[i,j] > qs[j]):
                condensationCount = condensationCount+1;
                #delta[i,j] = 10*vapor_n[i,j]/qs[j];
                delta[i,j] = 0.5;
                vapor_n[i,j] = qs[j];
            else:
                delta[i,j] = 0;

            if(vapor_n[i,j] < 0):
                vapor_n[i,j] = min(qs);
                negativeRHCount = negativeRHCount+1;
                ''' hack to ensure that RH is non-negative always.
                ''' 

    dt = tn - prev_t;

    ''' condensation strength is proportional to eta height'''
    #print 'delta: ', delta;
    #print 'height field: ', psi;
    delta = delta*etaDiff;
    #print 'forcing: ', delta;
    #raw_input();
    print 'forcing strength: ',sum(delta);


    (tn2,omega_n) = stepper_omega.integrate(tn,omega_n,(delta));
    print  'v,q coherence: ', sum(v*vapor);
    print 'Time elapsed: ',tn;
    print 'Number of cells condensing: ', condensationCount;
    print 'negative RH cells: ', negativeRHCount;
    print 'minimum vapor value: ', amin(amin(vapor));

    if(tn > sec):
        rh = zeros((N,N/2+1));
        for i in arange(N):
            for j in arange(N/2+1):
                rh[i,j] = (vapor_n[i,j]/qs[j]);

        (n,bins) = \
                   histogram(rh,bins=[-0.05,0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95,1.05],normed=True);
        handle.set_ydata(n);
       # handle.autoscale();
       # c1.update_bruteforce(handle);
        np.save('/home/joymm/VaporS'+str(beta)+'/rh_data_'+str(sec)+'.npy',rh.transpose());
        np.save('/home/joymm/VaporS'+str(beta)+'/delta_data_'+str(sec)+'.npy',delta.transpose());
        np.save('/home/joymm/VaporS'+str(beta)+'/vapor_data_'+str(sec)+'.npy',vapor_n.transpose());
        handle3.set_array(rh.transpose());
        handle3.autoscale();
        np.save('/home/joymm/VaporS'+str(beta)+'/eta_data_'+str(sec)+'.npy',omega_n.transpose());
        handle2.set_array(omega_n.transpose());
        handle2.autoscale();
        handle4.set_array(u.transpose());
        handle4.autoscale();
        handle5.set_array(v.transpose());
        handle5.autoscale();
        sec = sec+1;
        savefig('/home/joymm/VaporS'+str(beta)+'/fig'+str(sec)+'.png');
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
