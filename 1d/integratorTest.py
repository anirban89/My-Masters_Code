import numpy
from scipy import amax;


class integratorTest:

    dx = 0.0;
    callback = 0.0;

    def __init__(self,delta, advect, diffuse, maxvel, dim=1):

        """ Init routine for the integrator.
        
        delta: stores the spatial grid width

        advect: list which stores the functions to be used to calculate dm/dt
        for the advection step

        diffuse: list of functions which is called to do diffusion after the
        advection step

        maxvel: used to give the integrator information about the maximum
        velocity for calculating the time step

        dim: either 1 or 2 for now

        """

        print 'Initializing RK4 integrator';
        self.dx = delta;
        self.cflConstant = 0.5;
        if(dim == 2):
            print 'In 2D mode';
            self.dudt = advect[0];
            self.dvdt = advect[1];
            self.diffusionu = diffuse[0];
            self.diffusionv = diffuse[1];
            self.integrate = self.integrate2D;
            self.maxvelocity = maxvel;
        elif (dim == 1):
            print 'In 1D mode';
            self.dfdt = advect[0];
            self.diffusion = diffuse[0];
            self.integrate = self.integrate1D;
            self.maxvelocity = maxvel;
        print 'Done';

    def integrate1D(self,t,f,args=None):
        # integrate only one step in time.
        # assume same delta in x and y

        ''' f will normally be a stream function. The velocities are calculated
        in the callback routine, and the maxvelocity function, so the integrator
        only passes f along to these.
        '''

        maxf = self.maxvelocity(f, args);

        ''' get the time step from CFL condition and use half of it to run
        advection and the other half to do diffusion
        '''

        dT = self.cflConstant*self.dx/maxf;
        if(dT > 1):
            dT = 0.75;
        dt = dT/2;

        k1 = self.dfdt(t,f,args);

        k2 = self.dfdt(t+dt/2, f+(dt*k1/2.0), args);

        k3 = self.dfdt(t+dt/2, f+(dt*k2/2), args);

        k4 = self.dfdt(t+dt, f+(dt*k3), args);


        k = (k1 + 2*k2 + 2*k3 + k4)/6;

        fn = f + dt*k;

        ''' diffusion must be done for half the time step
        '''
        fn = self.diffusion(dt,fn,args);

        tn = t + dT;

        return (tn,fn);

       
    def integrate2D(self,t,u,v):
        # integrate only one step in time.
        # assume same delta in x and y
        maxVel = self.maxvelocity(u,v);

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


