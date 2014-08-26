from scipy import fftpack
from scipy import *

class diffs1d:

    dataLength = 0;

    def __init__(self, length):

        N = length;
        self.length = length;
        dataRange = linspace(0,N-1,N);
        self.h = 2*pi/length;

        print 'Initializing diffs\n';

        #DATA length must be even
        multiplier = 1.j*dataRange;

        self.firstOrder = (multiplier);
        self.secondOrder = (multiplier*multiplier);


    def partial (self, field, dim=0, order=1, noIfft=False):

        # order must be one or two, field must be 2 dimensional

        #by default, it calculates the first order partial along x axis

        fHat = fft(field);

        if order == 1:
            multiplier = self.firstOrder;
        elif order == 2:
            multiplier = self.secondOrder;

        print 'length: ',shape(multiplier);

        print 'fHat', shape(fHat), 'mult', shape(multiplier);
        fHat = fHat*multiplier;

        fPrime = ifft(fHat);

#return 2*real(fPrime)/(pow(self.length,order));
        return real(fPrime);

    def partialinv (self, field, dim=0, order=1, noIfft=False):

        # order must be one or two, field must be 2 dimensional

        #by default, it calculates the first order partial along x axis

        fHat = fft(field);

        if order == 1:
            multiplier = self.firstOrder;
        elif order == 2:
            multiplier = self.secondOrder;

        multiplier[0] = 1;

        print 'length: ',shape(multiplier);

        print 'fHat', shape(fHat), 'mult', shape(multiplier);
        fHat = fHat/multiplier;


        fPrime = ifft(fHat);

        return real(fPrime)*2;


    def partialSine(self, field, order=1):

        [lenX,lenY] = shape(field);
        zero = zeros((1,lenY));

        print 'zero', shape(zero), 'field', shape(field);
        tmp = \
        concatenate((zero,field,zero,-field[::-1]));
        print 'length: ',shape(tmp);
        
        result = self.partial(tmp,0);

        return result[1:len(field)+1];
