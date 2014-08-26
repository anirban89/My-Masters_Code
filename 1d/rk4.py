import numpy
from scipy import amax;


class rk4:

    dx = 0.0;
    callback = 0.0;

    def __init__(self,delta, dudt, dvdt):

        print 'Initializing RK4 integrator';
        self.dx = delta;
        self.cflConstant = 0.5;
        self.dudt = dudt;
        self.dvdt = dvdt;
        print 'Done';

    def integrate(self,t,u,v):
        # integrate only one step in time.
        # assume same delta in x and y
        maxu = amax(u);
        maxv = amax(v);
        maxVel = amax((maxu,maxv));

        dt = self.cflConstant*self.dx/maxVel;
        print 'Time step selected: ', dt;

        k1 = self.dudt(t,u,v);
        l1 = self.dvdt(t,u,v);

        k2 = self.dudt(t+dt/2, u+(dt*k1/2), v+(dt*l1/2));
        l2 = self.dvdt(t+dt/2, u+(dt*k1/2), v+(dt*l1/2));

        k3 = self.dudt(t+dt/2, u+(dt*k2/2), v+(dt*l2/2));
        l3 = self.dvdt(t+dt/2, u+(dt*k2/2), v+(dt*l2/2));

        k4 = self.dudt(t+dt, u+(dt*k3), v+(dt*l3));
        l4 = self.dvdt(t+dt, u+(dt*k3), v+(dt*l3));

        k = (k1 + 2*k2 + 2*k3 + k4)/6;
        l = (l1 + 2*l2 + 2*l3 + l4)/6;

        un = u + dt*k;
        vn = v + dt*k;
        tn = t + dt;

        return (tn,un,vn);



