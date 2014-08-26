import numpy
from scipy import amax,shape;


class integrator:

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

        dim: either 1 or 2  or 3 for now

        """

        print 'Initializing RK4 integrator';
        self.dx = delta;
        self.cflConstant = 0.5;
        if(dim == 3):
            print 'In 3D mode';
            self.dfdt = advect[0];
            self.diffusion = diffuse[0];
            self.integrate = self.integrate3D;
            self.maxvelocity = maxvel;

        elif(dim == 2):
            print 'In 2D mode';
            self.dfdt = advect[0];
            self.diffusion = diffuse[0];
            self.integrate = self.integrate2D;
            self.maxvelocity = maxvel;

        elif (dim == 1):
            print 'In 1D mode';
            self.dfdt = advect[0];
            self.diffusion = diffuse[0];
            self.integrate = self.integrate1D;
            self.maxvelocity = maxvel;
        print 'Done';

    def integrate1D(self,t,f,args=None, dt=0):
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

        if(dt > 0):
            dT = dt;
        else:
            dT = self.cflConstant*self.dx/maxf;
            dT = min(dT,1.0);

        if(dT > 1):
           dT = 1.0;

        dt = dT/2;

        print 'Time step selected for 1D: ', dT;

        k1 = self.dfdt(t,f,args);

        k2 = self.dfdt(t+dt/2, f+((dt/2)*k1), args);

        k3 = self.dfdt(t+dt/2, f+((dt/2)*k2), args);

        k4 = self.dfdt(t+dt, f+(dt*k3), args);

        k = (k1 + 2.0*k2 + 2.0*k3 + k4)/6;

        fn = f + dt*k;

        ''' diffusion must be done for half the time step
        '''
        fn = self.diffusion(dt,fn,args);


        return (dT,fn);
       
    def integrate2D(self,t,f,args=None, dt=0):

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

        if(dt > 0):
            dT = dt;
        else:
            dT = self.cflConstant*self.dx/maxf;
            dT = min(dT,0.1);

        if(dT > 1):
           dT = 1.0;

        dt = dT/2;

        print 'Time step selected for 1D: ', dT;

        k1 = self.dfdt(t,f,args);

        [x,y] = k1;
        k1n = [((dt/2)*x)+f[0],((dt/2)*y)+f[1]];

        k2 = self.dfdt(t+dt/2, k1n, args);

        [x,y] = k2;
        k2n = [((dt/2*x))+f[0],((dt/2)*y)+f[1]];

        k3 = self.dfdt(t+dt/2, k2n, args);

        [x,y] = k3;
        k3n = [(dt*x)+f[0],(dt*y)+f[1]];

        k4 = self.dfdt(t+dt, k3n, args);

        [x,y] = k1;
        k1 = [x/6,y/6];

        [x,y] = k2;
        k2 = [x/3,y/3];

        [x,y] = k3;
        k3 = [x/3,y/3];

        [x,y] = k4;
        k4 = [x/6,y/6];


        k = (k1 + k2 + k3 + k4);

        fn = [f[0]+dt*k[0], f[1]+dt*k[1]];

        ''' diffusion must be done for half the time step
        '''
        fn = self.diffusion(dt,fn,args);

        tn = t + dT;

        return (tn,fn);
 
       
    def integrate3D(self,t,f,args=None, dt=0):
        # integrate only one step in time.
        # assume same delta in x and y

        
        if(dt > 0):
            dT = dt;
        else:
            maxVel = self.maxvelocity((f), args);

            if(maxVel == 0):
                dT = 1.0;
            else:
                dT = self.cflConstant*self.dx/maxVel;

        dT = min((dT,1.0));
        dt = dT/2;

        print 'Time Step selected for 3D: ', dT;

        k1 = self.dfdt(t,f,args);

        [x,y,z] = k1;
        k1n = [((dt/2)*x)+f[0],((dt/2)*y)+f[1],((dt/2)*z)+f[2]];

        k2 = self.dfdt(t+dt/2, k1n, args);

        [x,y,z] = k2;
        k2n = [((dt/2)*x)+f[0],((dt/2)*y)+f[1],((dt/2)*z)+f[2]];

        k3 = self.dfdt(t+dt/2, k2n, args);

        [x,y,z] = k3;
        k3n = [(dt*x)+f[0],(dt*y)+f[1],(dt*z)+f[2]];

        k4 = self.dfdt(t+dt, k3n, args);

        [x,y,z] = k1;
        k1 = [x/6,y/6,z/6];

        [x,y,z] = k2;
        k2 = [x/3,y/3,z/3];

        [x,y,z] = k3;
        k3 = [x/3,y/3,z/3];

        [x,y,z] = k4;
        k4 = [x/6,y/6,z/6];


        k = (k1 + k2 + k3 + k4);

        fn = [f[0]+dt*k[0], f[1]+dt*k[1], f[2]+dt*k[2]];

        ''' diffusion must be done for half the time step
        '''
        fn = self.diffusion(dt,fn,args);


        return (dT,fn);


