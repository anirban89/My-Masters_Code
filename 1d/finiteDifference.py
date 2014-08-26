from pylab import *;

"""
This package implements finding partial derivatives on a finite
difference grid. A staggered grid is assumed, with a variable
'u' on the vertical, 'v' on the horizontal and 'n' in the centre
of each cell, i.e,

x -- v -- x
|         |
u    n    u
|         |   
x -- v -- x  

and so on, where 'x' marks the actual grid points.

since a partial derivate of u along the x direction will give values
of the derivative at the positions of 'n', the resulting values will
have to be average suitably to get the value at the desired location.

This is the reason why the variable 'gridPosition' is present as an
argument.
"""

def partialXfd(array, bc, variable='u', gridPosition='u', length=None):

    if(length == None):
        length = 2*pi;

    N,M = shape(array);

    dx = length/N;

    if(variable == 'u'):

        if(gridPosition == 'u'):

            temp = bc.applyBoundaryEW(array);

            N,M = shape(temp);

            diff = (temp[0:N-1,:] - temp[1:N,:])/dx;

            N,M = shape(diff);

            #Horizontally average
            diff = (diff[0:N-1,:] + diff[1:N,:])/2;

        elif(gridPosition == 'v'):

            temp = bc.applyBoundaryNS(array);
            N,M = shape(temp);

            diff = (temp[0:N-1,:] - temp[1:N,:])/dx;

            N,M = shape(diff);
            #Vertically average
            diff = (diff[:,0:M-1] + diff[:,1:M])/2;

        else: #gridPosition = 'n'

            diff = (array[0:N-1,:] - \
                        array[1:N,:])/dx;


        return diff;

    if(variable == 'v'):

        if(gridPosition == 'u'):

            temp = bc.applyBoundaryEW(array);

            diff = (temp[0:len(temp)-2,:] - temp[1:len(temp)-1,:])/dx;

            #Vertically average
            diff = (diff[:,0:len(diff)-2] + diff[:,1:len(diff)-1])/2;

        elif(gridPosition == 'v'):

            temp = bc.applyBoundaryEW(array);

            diff = (temp[0:len(temp)-2,:] - temp[1:len(temp)-1,:])/dx;

            #Horizontally average
            diff = (diff[0:len(diff)-2,:] + diff[1:len(diff)-1,:])/2;

        else: #gridPosition = 'n'

            temp = bc.applyBoundaryEW(array);

            diff = (temp[0:len(temp)-2,:] - temp[1,len(temp)-1,:])/dx;

            #Vertically average
            diff = (diff[:,0:len(diff)-2] + diff[:,1:len(diff)-1])/2;

            #Horizontally average
            diff = (diff[0:len(diff)-2,:] + diff[1:len(diff)-1,:])/2;

        return diff;

    if(variable == 'n'):

        if(gridPosition == 'u'):

            temp = bc.applyBoundaryEW(array);

            diff = (temp[0:len(temp)-2,:] - temp[1:len(temp)-1,:])/dx;

        elif(gridPosition == 'v'):


            temp = bc.applyCompleteBoundary(array, 'x');

            diff = (temp[:,0:len(temp)-2] - temp[:,1:len(temp)-1])/dx;

            #Horinzontally average
            diff = (diff[0:len(diff)-2,:] + diff[1:len(diff)-1,:])/2;

            #Vertically average
            diff = (diff[:,0:len(diff)-2] + diff[:,1:len(diff)-1])/2;


        else: # gridPosition = n

            temp = bc.applyBoundaryEW(array);

            diff = (temp[0:len(temp)-2,:] - temp[1:len(temp)-1,:])/dx;

            #Horizontally average
            diff = (diff[0:len(diff)-2,:] + diff[1:len(diff)-1,:])/2;

        return diff;



def boydPartialX(array, Length=None):

        if(Length == None):
                Length = 2*pi;


        if(array.ndim == 2):
            (xLen, yLen) = shape(array);

        if(array.ndim == 1):
            xLen, = shape(array);


        dx = Length/xLen;

        index = arange(0,xLen);

        #Periodic Boundaries, so wrap around

        arrayPlusOne = array[(index+1)%xLen,:];
        arrayPlusTwo = array[(index+2)%xLen,:];
        arrayPlusThree = array[(index+3)%xLen,:];
        arrayPlusFour = array[(index+4)%xLen,:];

        arrayMinusOne = array[(index-1),:];
        arrayMinusTwo = array[(index-2),:];
        arrayMinusThree = array[(index-3),:];
        arrayMinusFour = array[(index-4),:];

        diff = (4.0/5.0*(arrayPlusOne-arrayMinusOne) +
            1.0/5.0*(arrayMinusTwo-arrayPlusTwo)
            + 4.0/105.0*(arrayPlusThree-arrayMinusThree)
            + 1.0/280.0*(arrayMinusFour-arrayPlusFour))/dx;

        return diff;

def boydPartialY(array, bc, Length=None):

        if(Length == None):
                Length = 2*pi;


        if(array.ndim == 2):
            (xLen, yLen) = shape(array);
            temp = bc.applyBoundaryNS(array);

        if(array.ndim == 1):
            yLen, = shape(array);
            temp = bc.applyBoundaryNS1d(array);

        dy = Length/yLen;


        index = arange(0,yLen);

        arrayPlusOne = temp[:,(index+5)];
        arrayPlusTwo = temp[:,(index+6)];
        arrayPlusThree = temp[:,(index+7)];
        arrayPlusFour = temp[:,(index+8)];

        arrayMinusOne = temp[:,(index+3)];
        arrayMinusTwo = temp[:,(index+2)];
        arrayMinusThree = temp[:,(index+1)];
        arrayMinusFour = temp[:,(index)];

        diff = (4.0/5.0*(arrayPlusOne-arrayMinusOne) +
            1.0/5.0*(arrayMinusTwo-arrayPlusTwo)
            + 4.0/105.0*(arrayPlusThree-arrayMinusThree)
            + 1.0/280.0*(arrayMinusFour-arrayPlusFour))/dy;

        return diff;

