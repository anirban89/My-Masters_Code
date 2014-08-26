from coupledDerivative import *;
from integrator import *;
from pylab import *;
import time;


N = 512;
j = arange(0,N);
h = 2*pi/N;

inp = zeros((N,4));

for i in arange(4):
    inp[:,i] = sin(2*pi*j/N);

deriv = coupledDerivative(N,4,2*pi,2*pi);

def dudt(dt,u,args):

    return -u*deriv.partialX(u);

def diffusion(dt, u, args):
    return u;

def calcMax_omega(u,args):

    return amax(amax(u));

stepper_omega = integrator(h, [dudt], [diffusion], calcMax_omega);

ion();
output, =  plot(inp[:,1]);
#hold(True);
#plot(deriv.partialX(inp)[:,1]);
#draw();
#raw_input();
#exit();

while(1):

    print shape(inp);
    (dt, inp) = stepper_omega.integrate(0, inp,(),0);
    print shape(inp);
    output.set_ydata(inp[:,1]);
    draw();
    time.sleep(.1);
