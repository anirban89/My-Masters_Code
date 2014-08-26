from spectral import *;
from integrator import *;
from pylab import *;
from diffusion import *;
import scipy as sp;
import pylab as m;
import sys;
from os.path import expanduser;

home = expanduser("~");

resume = 0;
path = home+'/joymm/TwoLayer';

if(len(sys.argv) > 1):
    if(len(sys.argv) == 2):
    
        print "Resuming from";
        resume = int(sys.argv[1]);
        print resume;

N = 256;
h = 2*pi/N;

y = arange(0,N)-N/2;
y = y*h;

print y;

tn = 0;
sec = 1;

H = 3;
g = 10;
delta = 0.1;

u = 0.5*zeros((N,N));
v = 0.5*zeros((N,N));
h = rand(N,N);

def dudt(dt, u, args):
    """ 
        (psi_xx + psi_yy - psi) = omega. Invert, calculate u = -psi_y, v = psi_x
        omega_t = -(u*omega_x + v*omega_y + v + delta*psi)
    """
    (v,h) = args;
    vy = zeros((N,N));
    for i in arange(N):
        vy[i,:] = v[i,:]*y[i];

    return real(-(u*partialX(u) + v*partialY(u) - v - delta*vy + partialX(h)));

def dvdt(dt, v, args):
    """ 
        (psi_xx + psi_yy - psi) = omega. Invert, calculate u = -psi_y, v = psi_x
        omega_t = -(u*omega_x + v*omega_y + v + delta*psi)
    """
    (u,h) = args;
    uy = zeros((N,N));
    for i in arange(N):
        uy[i,:] = u[i,:]*y[i];

    return real(-(u*partialX(v) + v*partialY(v) + u + delta*uy + partialY(h)));


def calcMax(fields, args):
    (u,v,h) = fields;

    maximum = sqrt(amax(g*H,(amax(u*u),amax(v*v))));
    print 'Max vel: ', maximum;
    return maximum;

def dhdt(dt, h, args):
    """ 
        vapor_t = -(u*vapor_x + v*vapor_y + Beta*v)
    """
    (u,v,forcing) = args;

    hn = H+h;

    return real(-(partialX(hn*u) + partialY(hn*v) + forcing));


def dSystemdt(dt,fields,args):
    [u,v,h] = fields;

    (forcing) = args;

    dh = dhdt(dt, h, (u,v,forcing));
    du = dudt(dt, u, (v,h));
    dv = dvdt(dt, v, (u,h));
    return(du, dv, dh);


def diffusion(dt, fields, args):
    [u,v,h] = fields;

    diffh = spectralDiffusion(dt, h, args);
    diffu = spectralDiffusion(dt, u, args);
    diffv = spectralDiffusion(dt, v, args);

    return(du, dv, dh);


stepper = integrator(h, [dSystemdt], \
                [diffusion], calcMax_omega,3);

ion();
h_n = h;
u_n = u;
v_n = v;

dx = 2*pi/(N-1);
dy = 2*pi/(N-1);


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
   h_n = \
       np.load(path+'/h_data_'+str(resume)+'.npy').transpose();
   u_n = \
       np.load(path+'/u_data_'+str(resume)+'.npy').transpose();
   v_n = \
       np.load(path+'/v_data_'+str(resume)+'.npy').transpose();
   sec = resume+1;
   prev_t = resume;


subplot(121);
handle1 = imshow(v_n.transpose(), vmin=-15,vmax=1, aspect='auto');
#handle3 = imshow(vapor_n.transpose(), vmin=-15,vmax=1,cmap=my_cmap_vapor);
tick = linspace(-10,10,20);
c3 = colorbar();
subplot(122);
handle2 = imshow(u.transpose(), vmin=-15,vmax=1, aspect='auto');
c4 = colorbar();
draw();

forcingTime = 0;
global_time = 0;

while (sec<5000):

    forcing = 3.0*ones((N,N));
    #forcing = rand(N,N);
    (dt,fields) = stepper.integrate(tn,[u_n, v_n, \
                                            h_n],(forcing));

    [u_n, v_n,h_n] = fields;

    tn = dt + prev_t;
    prev_t = tn;
    print '============================='
    print 'Time Elapsed: ', tn;
    print '============================='


    if(tn > forcingTime):
        forcingTime = forcingTime+5;
        forcing = rand(N,N);

    if(tn > sec):

        print '==================================================================='

        print 'Time elapsed: ',tn;
        print 'negative RH cells: ', negativeRHCount;


        print '==================================================================='

        handle2.set_array(v_n.transpose());
        handle2.autoscale();
        handle1.set_array(u_n.transpose());
        handle1.autoscale();

        sec = sec+1;
        draw();
        np.save(path+'/v_data_'+str(sec)+'.npy',v_n.transpose());
        np.save(path+'/u_data_'+str(sec)+'.npy',u_n.transpose());
        np.save(path+'/h_data_'+str(sec)+'.npy',h_n.transpose());
        savefig(path+'/fig'+str(sec)+'.png');

ioff();

