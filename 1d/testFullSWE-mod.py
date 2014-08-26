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

y = zeros(N);

for i in arange(N):
    y[i] = (i-N/2+1)*h;


tn = 0;
sec = 1;
epsInv = 10;

'''max eta that should be reached'''
eta_u = 1;

#v = sin(t)+2*sin(t/4)+3*cos(t);
u = zeros((N,N));
v = zeros((N,N));
n = 0.001*ones((N,N));

#for i in arange(N):
#    for j in arange(N):
#        n[i,j] = 0.1*exp(-(pow(i-N/2,2)+pow(j-N/2,2))/100);


omega = (u,v,n);

vapor = zeros((N,N));
#for i in arange(20,80):
#    for j in arange(20,80):
#        omega[i,j] = 1000*rand();
#        omega[N-i,j] = value;
#        omega[i,N-j] = value;



def dudt(dt, omega, args):
    """ 
        (psi_xx + psi_yy - psi) = omega. Invert, calculate u = -psi_y, v = psi_x
        omega_t = -(u*omega_x + v*omega_y + v + delta*psi)
    """
    (u,v,n) = omega;
    (delta,q) = args;

    vy = zeros((N,N));

    for i in arange(N):
        vy[:,i] = v[:,i]*y[i];


    #u_n = real(-(u*partialX(u) + v*partialY(u) + partialX(n)*epsInv) + epsInv*v );
    u_n = real(-(u*partialX(u) + v*partialY(u) + partialX(n)*epsInv) + vy + epsInv*v );
#    print 'Advection: ', amax(amax(abs(real(-(u*partialX(u) +
#    v*partialY(u))) )))
#    print 'Pressure: ', amax(amax(abs(real(-( partialX(n)*epsInv)))))
#    print 'Beta effect: ',  amax(amax(abs((vy))))
#    print 'F effect: ', amax(amax(abs(epsInv*v )));
#    print 'UN: ',amax(amax(abs(u_n)));

    return(u_n);

def dvdt(dt, omega, args):
    """ 
        (psi_xx + psi_yy - psi) = omega. Invert, calculate u = -psi_y, v = psi_x
        omega_t = -(u*omega_x + v*omega_y + v + delta*psi)
    """
    (u,v,n) = omega;
    (delta,q) = args;

    uy = zeros((N,N));

    for i in arange(N):
        uy[:,i] = u[:,i]*y[i];

    #v_n = real(-(u*partialX(v) + v*partialY(v) + partialY(n)*epsInv + epsInv*u) );
    v_n = real(-(u*partialX(v) + v*partialY(v) + partialY(n)*epsInv + uy + epsInv*u) );
#    print 'VN: ',amax(amax(abs(v_n)));

    return(v_n);


def dndt(dt, omega, args):
    """ 
        (psi_xx + psi_yy - psi) = omega. Invert, calculate u = -psi_y, v = psi_x
        omega_t = -(u*omega_x + v*omega_y + v + delta*psi)
    """
    (u,v,n) = omega;
    (delta,q) = args;
    length,breadth = shape(n);


#    n_n = real(-(u*partialX(n) + v*partialY(n) + (partialX(u) + \
#    partialY(v))*(n + epsInv) - delta));
    n_n = real(-(u*partialX(n) + v*partialY(n) + (partialX(u) + \
    partialY(v))*(n + epsInv) - 0.01*rand(length, length)));
    #n_n = real(-(partialX(n*u) + partialY(n*v) + (partialX(u) + partialY(v))*epsInv) );

#    print 'nn: ', n_n;
    return(n_n);



def diffusion_SWE(dt,vals,args):

    (u,v,n) = vals;

    (u_n, nu) = diffusion_omega(dt, u, 0);
    
    (v_n, nu) = diffusion_omega(dt, v, nu);
    (n_n, nu) = diffusion_omega(dt, n, nu);

    return(u_n,v_n,n_n);


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

#    omegaHat[xVal/6:5*xVal/6,yVal/6:5*yVal/6] = 0;

    k = array(range(N), dtype=complex128) - N/2;
    k = ifftshift(k);
    [KX, KY] = meshgrid(k,k);
    mu = 0.0001;
    mask = zeros(N);
    for i in arange(N/100):
        mask[i] = 1;
        mask[N-i-1] = 1;

    [dragx,dragy] = meshgrid(mask,mask);

#   delsq = -(KX*KX + KY*KY)
    delsq = -(pow(KX,6) + pow(KY,6))
    omegasq = sqrt(sum(omega*omega)/(N*N));


    if (args == 0):
        nu = omegasq/pow(N,4);
    else:
        nu = args;
#    print 'viscosity=',nu;
    
    delsq = (nu*delsq - mu*(dragy + dragx))*dt;
    multiplier = exp(delsq);
    omegaHat = omegaHat*multiplier;
    return (real(ifft2(omegaHat)), nu);

def calcMax_omega(omega, args):
    (u,v,n) = (omega);
    g = 10;


    maximum = sqrt(amax((amax(u*u),amax(v*v))));
    print 'Max velocity: ', maximum;

    max_eta = amax(amax(abs(n)));
    gravity_speed = sqrt(g*max_eta);

    print 'Gravity wave velocity: ', gravity_speed; 
    return amax((maximum,gravity_speed));

def dudt_vapor(dt, vapor, args):
    """ 
        vapor_t = -(u*vapor_x + v*vapor_y + Beta*v)
    """
    (beta,u,v,delta) = args;

    return real(-(u*partialX(vapor) + v*partialY(vapor) - beta*v + delta*vapor));

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
    uMax = amax(abs(u));
    vMax = amax(abs(v));
    return amax((uMax,vMax));
#    maximum = sqrt(amax(amax(u*u),amax(v*v)));
#    return maximum;


def calcPower(field):
    [x,y] = shape(field);
    power = zeros(y);

    for i in arange(x):
        power = power + abs(fft(field[i,:]));

    power = power/x;
    return log(power);

stepper_omega = integrator(h, [dudt,dvdt,dndt], [diffusion_SWE], calcMax_omega, dim=3);
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
coeff = 0.1;

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


figure(figsize= (12,12));

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

(u,v,n) = omega;
prev_t = 0;
delta = zeros((N,N));
tn2 = 0;
sec = 0;
#jet();

if(resume > 0):
   print 'Resuming from step: ',resume;
   eta = \
   np.load('/home/joy/SWE'+str(coeff)+'/eta_data_'+str(resume)+'.npy').transpose();

   u = np.load('/home/joy/SWE'+str(coeff)+'/uVel_data_'+str(resume)+'.npy').transpose();

   v = \
   np.load('/home/joy/SWE'+str(coeff)+'/vVel_data_'+str(resume)+'.npy').transpose();

   vapor_n = \
   np.load('/home/joy/SWE'+str(coeff)+'/vapor_data_'+str(resume)+'.npy').transpose();

   omega_n = (u,v,eta);
   sec = resume+1;
   tn = resume;



subplot(221);
handle_u = imshow(u.transpose(),vmin=-0.1,vmax=0.1,cmap=my_cmap);
tick = linspace(-10,10,20);
c_u = colorbar();
draw();

subplot(222);
handle_v = imshow(v.transpose(),vmin=-0.1,vmax=0.1,cmap=my_cmap);
tick = linspace(-10,10,20);
c_v = colorbar();
draw();

subplot(223);
handle_n = imshow(n.transpose(),vmin=-0.1,vmax=0.1,cmap=my_cmap);
tick = linspace(-10,10,20);
c_n = colorbar();
draw();

subplot(224);
handle_vap = imshow(vapor_n.transpose(), vmin=-15,vmax=1);
#handle3 = imshow(vapor_n.transpose(), vmin=-15,vmax=1,cmap=my_cmap_vapor);
tick = linspace(-10,10,20);
c_vap = colorbar();
draw();


DT = 0.0005;

count = 0;

while (1):
    (u,v,n) = (omega_n);
    '''ensure that psi (eta) does not cross the ranges under which shallow water
    equations are valid'''
    count = count+1;

    #n[n>eta_u] = eta_u;
    #n[n<-eta_u] = -eta_u;


    condensationCount = 0;
    negativeRHCount = 0;
    vorticity = partialX(v) - partialY(u);
    (tn,vapor_n) = stepper_vapor.integrate(tn,vapor_n,(coeff,u,v,delta), dt=DT);

    '''actual value to multiply in momentum equation if condensation occurs'''
    etaDiff = eta_u - n;

    ''' if psi is greater than eta_u, then do not force system '''
    for i in arange(N):
        for j in arange(N):
            if(abs(n[i][j]) > eta_u):
                etaDiff[i][j] = 0;

    for i in arange(N):
        for j in arange(N):
            vapor[i,j] = vapor_n[i,j];
            delta[i,j] = 0;

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


    ''' condensation strength is proportional to eta height'''
    #print 'delta: ', delta;
    #print 'height field: ', psi;
    delta_eta = delta*etaDiff;
    #print 'forcing: ', delta;
    #raw_input();
#    print 'forcing strength: ',sum(delta);

    (DT,omega_n) = \
    stepper_omega.integrate(prev_t,omega_n,(delta_eta,vapor_n), 0);
    #DT = tn2 - prev_t;
    #tn = tn + DT;
    (u,v,n) = omega_n;
    print 'Time elapsed: ',tn;
    print 'Number of cells condensing: ', condensationCount;
#    print 'negative RH cells: ', negativeRHCount;
#    print 'minimum vapor value: ', amin(amin(vapor));

    if(count > 10):
       # handle.autoscale();
       # c1.update_bruteforce(handle);
       # np.save('/home/joy/Vapor'+str(beta)+'/rh_data_'+str(sec)+'.npy',rh.transpose());
       # np.save('/home/joy/Vapor'+str(beta)+'/vapor_data_'+str(sec)+'.npy',vapor_n.transpose());
       # np.save('/home/joy/Vapor'+str(beta)+'/delta_data_'+str(sec)+'.npy',delta.transpose());
        count = 0;
        handle_vap.set_array(vapor_n.transpose());
        handle_vap.autoscale();

        handle_u.set_array(u.transpose());
        handle_u.autoscale();

        handle_v.set_array(v.transpose());
        handle_v.autoscale();

        handle_n.set_array(n.transpose());
        handle_n.autoscale();
#        handle2.set_array(omega_n.transpose());
#        handle2.autoscale();
#        c2.update_bruteforce(handle2);
    if(tn > sec):
            savefig('/home/joy/SWE'+str(coeff)+'/fig'+str(sec)+'.png');
            np.save('/home/joy/SWE'+str(coeff)+'/eta_data_'+str(sec)+'.npy',n.transpose());
            np.save('/home/joy/SWE'+str(coeff)+'/uVel_data_'+str(sec)+'.npy',u.transpose());
            np.save('/home/joy/SWE'+str(coeff)+'/vVel_data_'+str(sec)+'.npy',v.transpose());
            np.save('/home/joy/SWE'+str(coeff)+'/vapor_data_'+str(sec)+'.npy',vapor_n.transpose());
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
