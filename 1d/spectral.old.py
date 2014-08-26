from numpy.fft import fft, ifft, ifftshift;
from scipy import *;

class spectral:

    dataLength = 0;

    def __init__(self, length=0):

        self.length = length;
        N = length;
        dataRange = linspace(0,N-1,N);

        print 'Initializing diffs\n';

        #DATA length must be even
        multiplier = 1.j*dataRange;

        self.firstOrder = (multiplier);
        self.secondOrder = (multiplier*multiplier);

        print 'Done Initializing diffs\n';


    '''
    def partial (self, field, dim=0, order=1, noIfft=False):

        # order must be one or two, field must be 2 dimensional
        #by default, it calculates the first order partial along x axis

        [xIndex,yIndex] = shape(field);

        fHat = fft(field, axis=dim);

        if order == 1:
            multiplier = self.firstOrder;
        elif order == 2:
            multiplier = self.secondOrder;

        print 'length: ',shape(multiplier);

        print 'fHat', shape(fHat), 'mult', shape(multiplier);
        if dim == 0:
           for i in arange(yIndex):
                fHat[:,i] = fHat[:,i]*(multiplier);

        elif dim == 1:
            for i in arange(xIndex):
                fHat[i,:] = fHat[i,:]*(multiplier);

        else:
            print 'invalid argument';
            return;

        if(noIfft):
            return fHat;

        fPrime = ifft(fHat,axis=dim);

        return real(fPrime);

    def partialSine(self, field, order=1):

        [lenX,lenY] = shape(field);
        zero = zeros((1,lenY));

        print 'zero', shape(zero), 'field', shape(field);
        tmp = \
        concatenate((zero,field,zero,-field[::-1]));
        print 'length: ',shape(tmp);
        
        result = self.partial(tmp,0);

        return result[1:len(field)+1];
    '''

    def InvPotentialVorticity(self,field, length=None):

        if length is None:
            length = 2*pi;

        N = shape(field)[0];

        k = array(range(N),dtype=complex128);
        k = concatenate((range(0,N/2),range(-N/2,0)));
        k *= (2*pi)/length;

        [KX, KY] = meshgrid(k,k);

        """ We are trying to solve d_yy(eta) + d_xx(eta) - eta = p
        Therefore, in Fourier space, it will become
        (-(kx^2 + ky^2) - 1)etaHat = pHat
        """
        delsq = -(KX*KX + KY*KY) - 1 ;
        delsq[0,0] = 1;

        tmp = fft(field,axis=0);
        tmp = fft(tmp,axis=1);

        tmp = tmp/delsq;

        tmp = ifft(tmp,axis=1);
        tmp = ifft(tmp,axis=0);

        return tmp.real;

    def InvLaplacian(self,field, length=None):

        if length is None:
            length = 2*pi;

        N = shape(field)[0];

        k = array(range(N),dtype=complex128);
        k = concatenate((range(0,N/2),range(-N/2,0)));
        k *= (2*pi)/length;

        [KX, KY] = meshgrid(k,k);

        """ We are trying to solve d_yy(eta) + d_xx(eta)  = p
        Therefore, in Fourier space, it will become
        (-(kx^2 + ky^2) )etaHat = pHat
        """
        delsq = -(KX*KX + KY*KY) ;
        delsq[0,0] = 1;

        tmp = fft(field,axis=0);
        tmp = fft(tmp,axis=1);

        tmp = tmp/delsq;

        tmp = ifft(tmp,axis=1);
        tmp = ifft(tmp,axis=0);

        return tmp.real;


    def partial(self, x, order=1, length=None, axis=0):
        """Numerically differentiate `x` using the pseudo-spectral method.

        Parameters
        ----------
        x : array_like
        The periodic data to be differentiated.
        order : int
        The order of the derivative to be computed.
        length : float
        The length of the domain on which the signal was sampled.
        axis : int
        The axis of `x` containing the data to be differentiated.

        Returns
        -------
        dx : array, with the same shape as `x`
        The differentiated data.
        """
        if length is None:
            length = 2 * pi;
        if axis is None:
            axis = -1;

        N = x.shape[axis];
        y = fft(x, axis=axis);

        k = array(range(N), dtype=complex128) - N/2;
        k *= 2*pi*1j / length;
        k = k**order;
        shp = ones(len(x.shape));
        shp[axis] = N;
        k.shape = shp;

        dy = ifftshift(k) * y;
        dx = ifft(dy, axis=axis).real;
        return dx;

    def partialX(field, order, length):
        return partial(field, order, length);
    
    def partialY(field, order, length):
        return partial(field, order, length, axis=1);

    def fourier(field,axis=0):
        return fft(field,axis);

    def invFourier(field,axis=0):
        return ifft(field,axis);
