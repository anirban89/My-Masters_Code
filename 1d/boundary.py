import numpy as np;
from scipy import concatenate, zeros,shape;


class boundary:
    
    east = 0;
    west = 0;
    north = 0; 
    south = 0;
    xboundaryType = 0;
    yboundaryType = 0;

    def __init__(self, xBcType='Periodic', yBcType='Periodic', bcs=None):
        """ Init for boundary adder

            the bcType can be 'Periodic', 'Neumann', 'AntiSymmetric', 'Symmetric' or 'Dirichlet'

            Cauchy boundary condition is not implemented.

            'Periodic' just adds the first grid as BC for the last and vice
            versa

            'Neumann' just extends the last grid point to one more point to
            ensure gradient at the end is zero.

            'Symmetric' mirrors four values at one end and 4 zeros at the
            other (From Boyd's implementation of Equatorial Beta
            plane)

            'AntiSymmetric' is same as above, except the values are
            negated to make it antisymmetric

            'Dirichlet' has specified values at the boundary, specified in the
            variable 'bcs'


        """

        self.xBoundaryType = xBcType;
        self.yBoundaryType = yBcType;

        if((self.xBoundaryType == 'Periodic')\
                 or (self.xBoundaryType == 'Neumann')):

            print 'No explicit zonal boundary conditions necessary'

        if((self.yBoundaryType == 'Periodic')\
                 or (self.yBoundaryType == 'Neumann')
                 or (self.yBoundaryType == 'AntiSymmetric')
                 or (self.yBoundaryType == 'Symmetric')):

            print 'No explicit meridional boundary conditions necessary'

        if(self.xBoundaryType == 'Dirichlet'):
            (self.east,self.west) = bcs;

        if(self.yBoundaryType == 'Dirichlet'):
            (self.north,self.south) = bcs;

        print 'Initialized boundary conditions';
        return;


    def applyBoundaryNS(self, array):

        (xLen, yLen) = shape(array);

        North = np.zeros((xLen,yLen));
        South = np.zeros((xLen,yLen));

        if(self.yBoundaryType == 'Neumann'):
           North[:,0] = array[:,yLen-1];
           South[:,0] = array[:,0];

        elif(self.yBoundaryType == 'Periodic'):
           South[:,0] = array[:,yLen-1];
           North[:,0] = array[:,0];

        elif(self.yBoundaryType == 'Symmetric'):
           South = np.zeros((xLen,4));
           North = np.zeros((xLen,4));

           South = array[:,1:5];
           #Flip it over, to make it symmetric
           South = South.transpose()[::-1].transpose();

        elif(self.yBoundaryType == 'AntiSymmetric'):
           South = np.zeros((xLen,4));
           North = np.zeros((xLen,4));

           South = array[:,1:5];
           #Flip it over, to make it symmetric
           South = -South.transpose()[::-1].transpose();


        else: # Dirichlet

           South = self.south;
           North = self.north;


        return(concatenate((South,array,North),axis=1));

    def applyBoundaryNS1d(self, array):

        yLen, = shape(array);

        North = np.zeros((yLen));
        South = np.zeros((yLen));

        South = self.south;
        North = self.north;


        return(concatenate((South,array,North)));


    def applyBoundaryEW(self, array):

        (xLen, yLen) = shape(array);

        East = np.zeros((1,yLen));
        West = np.zeros((1,yLen));

        if(self.xBoundaryType == 'Neumann'):
           East[0,:] = array[xLen-1,:];
           West[0,:] = array[0,:];

        elif(self.xBoundaryType == 'Periodic'):

           West[0,:] = array[xLen-1,:];
           East[0,:] = array[0,:];

        else:
           West[0,:] = self.west;
           East[0,:] = self.east;
            
        return(concatenate((West, array, East)));


    def returnBoundaries(self, array):

        (xLen, yLen) = shape(array);

        if(self.xBoundaryType == 'Neumann'):

           East = array[xLen-1,:];
           West = array[0,:];

        elif(self.xBoundaryType == 'Periodic'):

           West = array[xLen-1,:];
           East = array[0,:];

        else:
            West = self.west;
            East = self.east;

        if(self.yBoundaryType == 'Neumann'):

           North = array[:,yLen-1];
           South = array[:,0];

        elif(self.yBoundaryType == 'Periodic'):

           South = array[:,yLen-1];
           North = array[:,0];

        else:
           South = self.south;
           North = self.north;

        return(North, South, East, West);


    def applyCompleteBoundary(self, array, partialAlong= 'x'):

        '''
        Apply boundaries at all boundaries. Fill corners with
        appropriate boundary values.

        '''

        (xLen, yLen) = shape(array);

        (north, south, east, west) = self.returnBoundaries(array);

        temp = zeros((xLen+2, yLen+2));

        (xLen, yLen) = shape(temp);

        temp[1:xLen-1,1:yLen-1] = array[:,:];

        temp[0,1:xLen-1] = west[:];
        temp[xLen-1, 1:yLen-1] = east[:];

        temp[1:xLen-1,0] = south[:];
        temp[1:xLen-1,yLen-1] = north[:];


        if(partialAlong == 'x'):
            '''

                The north-south boundaries don't affect the corner
                boundaries. extend east-west boundaries to the
                corners.

            '''

            if(self.xBoundaryType == 'Periodic'):

                 temp[0,0] = temp[0,yLen-2];
                 temp[0,yLen-1] = temp[0,1];

                 temp[xLen-1,0] = temp[xLen-1,yLen-2];
                 temp[xLen-1,yLen-1] = temp[xLen-1,1];

            elif(self.xBoundaryType == 'Neumann'):

                 temp[0,0] = temp[0,1];
                 temp[0,yLen-1] = temp[0,yLen-2];

                 temp[xLen-1,0] = temp[xLen-1,1];
                 temp[xLen-1,yLen-1] = temp[xLen-1,yLen-2];

            else:

                 temp[0,0] = temp[1,0];
                 temp[0,yLen-1] = temp[1,yLen-1];

                 temp[xLen-1,0] = temp[xLen-2,0];
                 temp[xLen-1,yLen-1] = temp[xLen-2,yLen-1];

        if(partialAlong == 'y'):
             '''
             The exact opposite applies

             '''

             if(self.yBoundaryType == 'Periodic'):

                 temp[0,0] = temp[xLen-2,0];
                 temp[0,xLen-1] = temp[1,0];

                 temp[xLen-1,0] = temp[xLen-2,yLen-1];
                 temp[xLen-1,yLen-1] = temp[1,yLen-1];

             elif(self.yBoundaryType == 'Neumann'):

                 temp[0,0] = temp[1,0];
                 temp[0,yLen-1] = temp[xLen-1,0];

                 temp[xLen-1,0] = temp[1,yLen-1];
                 temp[xLen-1,yLen-1] = temp[xLen-1,yLen-2];

             else:

                 temp[0,0] = temp[1,0];
                 temp[0,yLen-1] = temp[1,yLen-1];

                 temp[xLen-1,0] = temp[xLen-2,0];
                 temp[xLen-1,yLen-1] = temp[xLen-2,yLen-1];

        return temp;
