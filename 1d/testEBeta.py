from finiteDifference import *;
from integrator import *;
from boundary import *;
from pylab import *;
import scipy as sp;
import pylab as m;
import sys;

resume = 0;
path = '/home/joy/NewEquatorialBeta/';

if(len(sys.argv) > 1):
     print "Resuming";
     resume = int(sys.argv[1]);


N = 128;

g = 10;
H = 250;
beta = 10**(-10);
Ly = 8*sqrt(sqrt(g*H)/beta);
Lx = 4*(10**7);

hy = Ly/N;
hx = Lx/(2*N);


y = zeros(N);
y1 = zeros(N);
epsilon = zeros(N);

for i in arange(N):
    y1[i] = (i-N/2+1)*hy;
    if(i < N/4 or i > 3*N/4):
        epsilon[i] = .001;


u = zeros((2*N,N));
v = zeros((2*N,N));
n = zeros((2*N,N));

x = linspace(0,Lx,2*N)
y = linspace(0,Ly,N)
Q = zeros((2*N,N));

for j in arange(N):
    for i in arange(2*N):
        #Q[i][j] = -0.1*exp(-(pow(y[j]-4*pi,2)+pow(x[i]-4*pi,2))/(4*pi));
        #Q[i][j] = \
        #        -3*(10**-5)*exp(-(pow(y[j]-(Ly/2),2)+pow(x[i]-(Lx/4),2))/(2000*Lx));
        Q[i][j] = \
                -3*(10**-5)*exp(-(pow(y[j]-(Ly/2),2)+pow(x[i]-(Lx/4),2))/(2000*Lx));

#n = Q;
#u = -Q;
#imshow(Q);
#colorbar();
#show();

#for i in arange(N):
#    for j in arange(N):
#        n[i,j] = 0.1*exp(-(pow(i-N/2,2)+pow(j-N/2,2))/100);


omega = (u,v,n);

#for i in arange(20,80):
#    for j in arange(20,80):
#        omega[i,j] = 1000*rand();
#        omega[N-i,j] = value;
#        omega[i,N-j] = value;

d = 0.0005;

def dudt(dt, omega, args):
    """ 
        (psi_xx + psi_yy - psi) = omega. Invert, calculate u = -psi_y, v = psi_x
        omega_t = -(u*omega_x + v*omega_y + v + delta*psi)
    """
    (u,v,n) = omega;
    (bc, Qn) = args;

    vy = zeros((2*N,N));

    for i in arange(N):
        vy[:,i] = v[:,i]*y1[i];

    pressure = g*boydPartialX(n,Lx);
    print 'Rotation ', average(abs(beta*vy));
#    print 'V ', average(abs(v));
    print 'Pressure ',average(abs(pressure));


    #u_n = real(-(boydPartialX(n,16*pi) + epsilon*u) + vy/(4*pi));
#    u_n = real(-( v*boydPartialY(u,bc,8*pi) + d*u/epsilon) + vy);
    #u_n = real(-( v*boydPartialY(u,bc,8*pi) - d*u/epsilon) + vy);
    u_n = real(-(pressure - beta*vy ));

    return(u_n);

def dvdt(dt, omega, args):
    """ 
        (psi_xx + psi_yy - psi) = omega. Invert, calculate u = -psi_y, v = psi_x
        omega_t = -(u*omega_x + v*omega_y + v + delta*psi)
    """
    (u,v,n) = omega;
    (bc, Qn) = args;

    uy = zeros((2*N,N));
    epsV = zeros((2*N,N));

    for i in arange(N):
        uy[:,i] = u[:,i]*y1[i];
        epsV[:,i] = v[:,i]*epsilon[i];


    v_n = real(-(g*boydPartialY(n,bc,Ly) + beta*uy + d*v));
    #v_n = real(-(g*boydPartialY(n,bc,Ly) + beta*uy + d*v));

    return(v_n);


def dndt(dt, omega, args):
    """ 
        (psi_xx + psi_yy - psi) = omega. Invert, calculate u = -psi_y, v = psi_x
        omega_t = -(u*omega_x + v*omega_y + v + delta*psi)
    """
    (u,v,n) = omega;
    (bc, Qn) = args;


    epsN = zeros((2*N,N));
    for i in arange(N):
        epsN[:,i] = n[:,i]*epsilon[i];

    divergence = -H*(boydPartialY(v,bc,Ly) + boydPartialX(u,Lx));
    #n_n = real(-(boydPartialX(u,16*pi) + boydPartialY(v, bc, 8*pi) + epsilon*n + Qn));
#    n_n = real(-(boydPartialY(h*v,bc, 8*pi) + boydPartialY(v,bc,8*pi) + d*n/epsilon) + Q/epsilon);
    #n_n = real(-(boydPartialY(h*v,bc, 8*pi) + boydPartialY(v,bc,8*pi) + d*n) + Q);
    n_n = real( divergence + Qn);


    return(n_n);



def dSystemdt(dt, omega, args):

    u_n = dudt(dt, omega, args);
    v_n = dvdt(dt, omega, args);
    n_n = dndt(dt, omega, args);

    return [u_n, v_n, n_n];

def diffusion_SWE(dt,vals,args):

    return vals;



def calcMax_omega(omega, args):
    (u,v,n) = (omega);


    maximum = sqrt(amax((amax(u*u),amax(v*v))));
    print 'Max velocity: ', maximum;

    max_eta = amax(amax(abs(n)));
    print 'Max eta: ', max_eta;
    gravity_speed = sqrt(g*H);

    print 'Gravity wave velocity: ', gravity_speed; 
    return 4*amax((maximum,gravity_speed));


def calcPower(field):
    [x,y] = shape(field);
    power = zeros(y);

    for i in arange(x):
        power = power + abs(fft(field[i,:]));

    power = power/x;
    return log(power);

stepper_omega = integrator(hx, [dSystemdt], [diffusion_SWE], calcMax_omega, dim=3);

ion();


prev_t = 0;

omega_n = (u,v,n);
#jet();

if(resume > 0):
   print 'Resuming from step: ',resume;
   eta = \
   np.load(path+'/eta_data_'+str(resume)+'.npy');

   u = np.load(path+'/uVel_data_'+str(resume)+'.npy');

   v = \
   np.load(path+'/vVel_data_'+str(resume)+'.npy');

   omega_n = (u,v,eta);
   sec = resume+1;
   tn = resume;



ion();

figure(figsize=(18,6));
handle_n = imshow(n.transpose());
colorbar();
draw();


DT=0;

north = zeros((2*N,4));
south = zeros((2*N,4));
bc = boundary(yBcType='Dirichlet',bcs=(north,south));
numerator = zeros((2*N,N));

tn = 0;
nextUpdate = tn + 3600;
hour = 0;

while (True):
    (u,v,n) = (omega_n);
    '''ensure that psi (eta) does not cross the ranges under which shallow water
    equations are valid'''

    #n[n>eta_u] = eta_u;
    #n[n<-eta_u] = -eta_u;


#    print 'forcing strength: ',sum(delta);

    (DT,omega_n) = \
    stepper_omega.integrate(prev_t,omega_n,(bc, Q));
    prev_t = tn;
    print DT;
    tn = tn + DT/2.0;
    #DT = tn2 - prev_t;
    #tn = tn + DT;
    (u,v,n) = omega_n;
    print 'Time elapsed: ',tn/3600.0, 'hours';
#    print 'negative RH cells: ', negativeRHCount;
#    print 'minimum vapor value: ', amin(amin(vapor));
#    if(tn > 1):
#        Q = zeros((N));

    if(tn >= nextUpdate):


        nextUpdate = tn + 3600;
        hour = hour+1;

        handle_n.set_array(n.transpose());
        handle_n.autoscale();

        print '--------------------------------------'
        print 'Energy:', sum(u*u + v*v + 10*n*n);
        print '--------------------------------------'

#        handle2.set_array(omega_n.transpose());
#        handle2.autoscale();
#        c2.update_bruteforce(handle2);
        savefig(path+'/fig'+str(hour)+'.png');
        np.save(path+'/eta_data_'+str(hour)+'.npy',n);
        np.save(path+'/uVel_data_'+str(hour)+'.npy',u);
        np.save(path+'/vVel_data_'+str(hour)+'.npy',v);
        draw();

ioff();

