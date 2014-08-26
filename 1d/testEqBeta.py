from coupledDerivative import *;
from integrator import *;
from boundary import *;
from pylab import *;
import scipy as sp;
import pylab as m;
import sys;

resume = 0;
path = '/home/joy/EquBeta/';

if(len(sys.argv) > 1):
     print "Resuming";
     resume = int(sys.argv[1]);


N = 256;
t = linspace(0,2*pi,N+1);
t = t[0:len(t)-1];
h = 8*pi/N;

y = zeros(N);

for i in arange(N):
    y[i] = (i-N/2+1)*h;


tn = 0;
sec = 1;
epsInv = 10;

'''max eta that should be reached'''
eta_u = 1;

#v = sin(t)+2*sin(t/4)+3*cos(t);
u = zeros((2*N,N));
v = zeros((2*N,N));
n = 0.01*ones((2*N,N));
ind = linspace(0,N,N);
vals = exp((-pow(ind-N/2,2)/20));
x = linspace(0,16*pi,2*N)
y = linspace(0,8*pi,N)


for i in arange(2*N):
    for j in arange(N):
        n[i,j] = rand();
#        n[i,j] = \
#        -exp((-pow(x[i]-8*pi,2)/10))*exp(-pow(y[j]-4*pi,2)/10)*cos(i);


#plot(n);
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



def dudt(dt, omega, args):
    """ 
        (psi_xx + psi_yy - psi) = omega. Invert, calculate u = -psi_y, v = psi_x
        omega_t = -(u*omega_x + v*omega_y + v + delta*psi)
    """
    (u,v,n) = omega;
    (bc, diff) = args;

    vy = zeros((2*N,N));

    for i in arange(N):
        vy[:,i] = v[:,i]*y[i];


    u_n = real(-(u*diff.partialX(u) + v*diff.partialY(u,bc) + diff.partialX(n)*epsInv) + epsInv*vy);
#    u_n = real(-u*diff.partialX(u) - v*diff.partialY(u,bc) - 10*diff.partialX(n) + 0.001*vy);

    return(u_n);

def dvdt(dt, omega, args):
    """ 
        (psi_xx + psi_yy - psi) = omega. Invert, calculate u = -psi_y, v = psi_x
        omega_t = -(u*omega_x + v*omega_y + v + delta*psi)
    """
    (u,v,n) = omega;
    (bc, diff) = args;

    uy = zeros((2*N,N));

    for i in arange(N):
        uy[:,i] = u[:,i]*y[i];

    #v_n = real(-(u*diff.partialX(v) + v*diff.partialY(v) + diff.partialY(n)*epsInv + epsInv*u) );
    v_n = real(-(u*diff.partialX(v) + v*diff.partialY(v, bc) + diff.partialY(n, bc)*epsInv + epsInv*uy) );
#    v_n = real(-u*diff.partialX(v) - v*diff.partialY(v, bc) - 10*diff.partialY(n, bc) - 0.001*uy);
#    print 'VN: ',amax(amax(abs(v_n)));

    return(v_n);


def dndt(dt, omega, args):
    """ 
        (psi_xx + psi_yy - psi) = omega. Invert, calculate u = -psi_y, v = psi_x
        omega_t = -(u*omega_x + v*omega_y + v + delta*psi)
    """
    (u,v,n) = omega;
    (bc, diff) = args;


#    n_n = real(-(diff.partialX(u) + diff.partialY(v, bc)));
    n_n = real(-(diff.partialX(n*u) + diff.partialY(n*v, bc) + (diff.partialX(u) + diff.partialY(v, bc))*epsInv));
#    n_n = real(-(diff.partialX(n*u) + diff.partialY(n*v,bc) + (diff.partialX(u) + diff.partialY(v,bc))) );

#    print 'nn: ', n_n;
    return(n_n);



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

stepper_omega = integrator(h, [dudt,dvdt,dndt], [diffusion_SWE], calcMax_omega, dim=3);

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
handle_u = imshow(u.transpose(),vmin=-0.1,vmax=0.1);
tick = linspace(-10,10,20);
c_u = colorbar();
draw();

subplot(132);
handle_v = imshow(v.transpose(),vmin=-0.1,vmax=0.1);
tick = linspace(-10,10,20);
c_v = colorbar();
draw();

subplot(133);
handle_n = imshow(n.transpose(),vmin=-0.1,vmax=0.1);
tick = linspace(-10,10,20);
c_n = colorbar();
draw();

dT = 0.0001;
DT=0;

count = 0;
north = zeros((2*N,10));
south = zeros((2*N,10));
bc = boundary(yBcType='Dirichlet',bcs=(north,south));
numerator = zeros((2*N,N));

print 'Initializing Derivatives'
deriv = coupledDerivative(2*N,N,16*pi,8*pi);

while (1):
    (u,v,n) = (omega_n);
    '''ensure that psi (eta) does not cross the ranges under which shallow water
    equations are valid'''
    count = count+1;

    #n[n>eta_u] = eta_u;
    #n[n<-eta_u] = -eta_u;


#    print 'forcing strength: ',sum(delta);

    (DT,omega_n) = \
    stepper_omega.integrate(prev_t,omega_n,(bc, deriv), dT);
    prev_t = tn;
    tn = tn + DT;
    #DT = tn2 - prev_t;
    #tn = tn + DT;
    (u,v,n) = omega_n;
    print 'Time elapsed: ',tn;
#    print 'negative RH cells: ', negativeRHCount;
#    print 'minimum vapor value: ', amin(amin(vapor));

    if(count >= 1):
       # handle.autoscale();
       # c1.update_bruteforce(handle);
       # np.save('/home/joy/Vapor'+str(beta)+'/rh_data_'+str(sec)+'.npy',rh.transpose());
       # np.save('/home/joy/Vapor'+str(beta)+'/vapor_data_'+str(sec)+'.npy',vapor_n.transpose());
       # np.save('/home/joy/Vapor'+str(beta)+'/delta_data_'+str(sec)+'.npy',delta.transpose());
        count = 0;

        handle_u.set_array(u.transpose());
        handle_u.autoscale();

        handle_v.set_array(v.transpose());
        handle_v.autoscale();

        handle_n.set_array(n.transpose());
        handle_n.autoscale();


        print '--------------------------------------'
        print 'Energy:', sum(u*u + v*v + 10*n*n);
        omega = (deriv.partialX(v) - deriv.partialY(u,bc));
        for i in arange(N):
                numerator[:,i] = epsInv*y[i] + omega[:,i];
        denominator = epsInv + n;
        print 'PV:', sum(numerator/denominator);
        print 'Enstrophy:', sum(omega*omega);
        print '--------------------------------------'

#        handle2.set_array(omega_n.transpose());
#        handle2.autoscale();
#        c2.update_bruteforce(handle2);
    if(tn > sec):
            savefig('/home/joy/EqBeta/fig'+str(sec)+'.png');
            np.save('/home/joy/EqBeta/eta_data_'+str(sec)+'.npy',n.transpose());
            np.save('/home/joy/EqBeta/uVel_data_'+str(sec)+'.npy',u.transpose());
            np.save('/home/joy/EqBeta/vVel_data_'+str(sec)+'.npy',v.transpose());
            sec = sec+1;
    draw();

ioff();

