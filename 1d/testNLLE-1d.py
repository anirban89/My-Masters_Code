from finiteDifference import *;
from integrator import *;
from boundary import *;
from pylab import *;
import scipy as sp;
import pylab as m;
import sys;

resume = 0;
path = '/home/joymm/NLLE-1d/';

if(len(sys.argv) > 1):
     print "Resuming";
     resume = int(sys.argv[1]);


N = 1024;
t = linspace(0,2*pi,N+1);
t = t[0:len(t)-1];
h = 8*pi/N;

y = zeros(N);
y1 = zeros(N);

for i in arange(N):
    y1[i] = (i-N/2+1)*h;



tn = 0;
sec = 1;
epsInv = 10;
epsilon = 0.1;

'''max eta that should be reached'''
eta_u = 1;

#v = sin(t)+2*sin(t/4)+3*cos(t);
u = zeros((N));
v = zeros((N));
n = 0.0*ones((N));
ind = linspace(0,N,N);
vals = exp((-pow(ind-N/2,2)/20));
x = linspace(0,16*pi,2*N)
y = linspace(0,8*pi,N)
Q = zeros((N));

for j in arange(N):
        Q[j] = -0.1*exp(-pow(y[j]-4*pi,2)/(4*pi));

print y1
n = Q;
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

d = 0.005;

def dudt(dt, omega, args):
    """ 
        (psi_xx + psi_yy - psi) = omega. Invert, calculate u = -psi_y, v = psi_x
        omega_t = -(u*omega_x + v*omega_y + v + delta*psi)
    """
    (u,v,n) = omega;
    (bc, Qn) = args;

    vy = zeros((N));

    for i in arange(N):
        vy[i] = v[i]*y1[i];


    #u_n = real(-(boydPartialX(n,16*pi) + epsilon*u) + vy/(4*pi));
#    u_n = real(-( v*boydPartialY(u,bc,8*pi) + d*u/epsilon) + vy);
    #u_n = real(-( v*boydPartialY(u,bc,8*pi) - d*u/epsilon) + vy);
    u_n = real(-( v*boydPartialY(u,bc,8*pi)) + vy);

    print 'UN:', shape(u_n);
    return(u_n);

def dvdt(dt, omega, args):
    """ 
        (psi_xx + psi_yy - psi) = omega. Invert, calculate u = -psi_y, v = psi_x
        omega_t = -(u*omega_x + v*omega_y + v + delta*psi)
    """
    (u,v,n) = omega;
    (bc, Qn) = args;

    uy = zeros((N));

    for i in arange(N):
        uy[i] = u[i]*y1[i];

    epsSq = epsilon*epsilon;

    #v_n = real(-(boydPartialY(n, bc, 8*pi) + epsilon*v) - uy/(4*pi));
#    v_n = real(-(v*boydPartialY(v,bc,8*pi)) - (boydPartialY(n,bc,8*pi) + uy)/epsSq - d*u/epsSq);
    #v_n = real(-(v*boydPartialY(v,bc,8*pi)) - (boydPartialY(n,bc,8*pi) + uy + d*v)/epsSq);
    v_n = real(-(v*boydPartialY(v,bc,8*pi)) - (boydPartialY(n,bc,8*pi) + uy)/epsSq);

    return(v_n);


def dndt(dt, omega, args):
    """ 
        (psi_xx + psi_yy - psi) = omega. Invert, calculate u = -psi_y, v = psi_x
        omega_t = -(u*omega_x + v*omega_y + v + delta*psi)
    """
    (u,v,n) = omega;
    (bc, Qn) = args;


    #n_n = real(-(boydPartialX(u,16*pi) + boydPartialY(v, bc, 8*pi) + epsilon*n + Qn));
#    n_n = real(-(boydPartialY(h*v,bc, 8*pi) + boydPartialY(v,bc,8*pi) + d*n/epsilon) + Q/epsilon);
    n_n = real(-(boydPartialY(n*v,bc, 8*pi) + boydPartialY(v,bc,8*pi) + d*n/epsilon) + (Q-n));
    n_n = real(-(boydPartialY(n*v,bc, 8*pi) + boydPartialY(v,bc,8*pi) + d*n/epsilon));

    return(n_n);

def dOmegadt(dt, omega, args):

    unew = dudt(dt,omega,args);
    vnew = dvdt(dt,omega,args);
    nnew = dndt(dt,omega,args);

    return((unew,vnew,nnew));

def diffusion_SWE(dt,vals,args):

    return vals;

    return(u_n,v_n,n_n);



def calcMax_omega(omega, args):
    (u,v,n) = (omega);
    g = 10;


    maximum = sqrt(amax((amax(u*u),amax(v*v))));
    print 'Max velocity: ', maximum;

    max_eta = amax(amax(abs(n)));
    gravity_speed = sqrt(g*max_eta);

    print 'Gravity wave velocity: ', gravity_speed; 
    return amax((maximum,gravity_speed));


def calcPower(field):
    [x,y] = shape(field);
    power = zeros(y);

    for i in arange(x):
        power = power + abs(fft(field[i,:]));

    power = power/x;
    return log(power);

stepper_omega = integrator(h, [dOmegadt], [diffusion_SWE], calcMax_omega, dim=3);

ion();

figure(figsize= (18,6));

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

omega_n = (u,v,n);
#jet();

if(resume > 0):
   print 'Resuming from step: ',resume;
   eta = \
   np.load('/home/joy/EqBeta/eta_data_'+str(resume)+'.npy').transpose();

   u = np.load('/home/joy/EqBeta/uVel_data_'+str(resume)+'.npy').transpose();

   v = \
   np.load('/home/joy/EqBeta/vVel_data_'+str(resume)+'.npy').transpose();

   omega_n = (u,v,eta);
   sec = resume+1;
   tn = resume;



subplot(131);
handle_u, = plot(u);
ylim((-0.05,0.05)); 
draw();

subplot(132);
handle_v, = plot(v);
ylim((-0.5,.5)); 
draw();

subplot(133);
handle_n, = plot(n);
ylim((-0.5,0.5)); 
draw();

DT = 0.0005;

count = 0;
north = zeros((4));
south = zeros((4));
bc = boundary(yBcType='Dirichlet',bcs=(north,south));
print bc;
numerator = zeros((N));

while (tn < 2000):
    (u,v,n) = (omega_n);
    '''ensure that psi (eta) does not cross the ranges under which shallow water
    equations are valid'''
    count = count+1;

    #n[n>eta_u] = eta_u;
    #n[n<-eta_u] = -eta_u;


#    print 'forcing strength: ',sum(delta);

    (DT,omega_n) = \
    stepper_omega.integrate(prev_t,omega_n,(bc, Q), 0);
    prev_t = tn;
    tn = tn + DT;
    #DT = tn2 - prev_t;
    #tn = tn + DT;
    (u,v,n) = omega_n;
    print 'Time elapsed: ',tn;
#    print 'negative RH cells: ', negativeRHCount;
#    print 'minimum vapor value: ', amin(amin(vapor));
#    if(tn > 1):
#        Q = zeros((N));

    if(count >= 1):
       # handle.autoscale();
       # c1.update_bruteforce(handle);
       # np.save('/home/joy/Vapor'+str(beta)+'/rh_data_'+str(sec)+'.npy',rh.transpose());
       # np.save('/home/joy/Vapor'+str(beta)+'/vapor_data_'+str(sec)+'.npy',vapor_n.transpose());
       # np.save('/home/joy/Vapor'+str(beta)+'/delta_data_'+str(sec)+'.npy',delta.transpose());
        count = 0;

        handle_u.set_ydata(u);
        handle_v.set_ydata(v);
        handle_n.set_ydata(n);


        print '--------------------------------------'
        print 'Energy:', sum(u*u + v*v + 10*n*n);
        print '--------------------------------------'

#        handle2.set_array(omega_n.transpose());
#        handle2.autoscale();
#        c2.update_bruteforce(handle2);
    if(tn > sec):
            savefig(path+'/fig'+str(sec)+'.png');
            np.save(path+'/eta_data_'+str(sec)+'.npy',n);
            np.save(path+'/uVel_data_'+str(sec)+'.npy',u);
            np.save(path+'/vVel_data_'+str(sec)+'.npy',v);
            sec = sec+1;
    draw();

ioff();

