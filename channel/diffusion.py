from numpy.fft import fft, ifft, ifftshift,fft2,ifft2;
from scipy import *

'''
Integrate the diffusion equation, in spectral domain for now.

'''

def specDiffusion1D(dt, u):
    '''
        du/dt = u_xx 
        => duHat/uHat = -k^2*dt
        => uHat = e^(-k^2)dt*uHat_0
    '''
    uHat = fft(u);
    numvals = len(u);
    threshold = int(numvals/3);
    uHat[threshold:2*numvals/3] = 0;
    k = array(range(N), dtype=complex128) - N/2;
    k = ifftshift(k);
    k = diffCoeff*(k**2);
    k = -k*dt;
    multiplier = exp(k);
    uHat = uHat*multiplier;
    return real(ifft(uHat));

def spectralDiffusion(dt, omega, args):
    '''
        d(omega)/dt = nu*((omega)_xx + (omega)_yy) - mu*omega
        => d(omega)Hat/(omega)Hat = -nu*(k_x^2 + k_y^2)*dt - mu*dt
        => (omega)Hat = e^(-nu*|k|^2 - mu)dt*(omega)Hat_0

        Dynamic viscosity nu = 0.01 for now
    '''
    omegaHat = fft2(omega);
    [xVal,yVal] = shape(omegaHat);
    numvals = len(omega);
    N = len(omega);
    downfilter = exp((arange(numvals/6+1, (numvals/2)+1) -
    numvals/2.0))[::-1];
    upfilter = exp(-(arange(numvals/2, (5*numvals/6)+1) -
    numvals/2.0))[::-1];

    preserveVals = ones(numvals/6);

    fftFilter = concatenate((preserveVals, downfilter, upfilter, preserveVals))

    print len(fftFilter), 'in diffusion'

    omegaHat = (fftFilter * omegaHat.transpose()).transpose();
    omegaHat = (fftFilter * omegaHat);



    k = array(range(N), dtype=complex128) - N/2;
    k = ifftshift(k);
    [KX, KY] = meshgrid(k,k);

    mu = 0.0000;
    mask = zeros(N);
    for i in arange(int(N/100)):
        mask[i] = 1;
        mask[N-i-1] = 1;

    [dragx,dragy] = meshgrid(mask,mask);


    delsq = -(pow(KX,6) + pow(KY,6))
    omegasq = sqrt(sum(omega*omega)/(N*N));

    nu = 10*omegasq/pow(N,4);
    print 'viscosity=',nu;
    
    delsq = (nu*delsq - mu*(dragx+dragy))*dt;
    multiplier = exp(delsq);
    omegaHat = omegaHat*multiplier;
    return real(ifft2(omegaHat));

