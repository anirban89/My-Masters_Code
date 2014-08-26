import spectral
import spectralfftw
from scipy import rand
from numpy.fft import fft, ifft
import anfft

a = rand(200,200);


scipy_a = spectral.partialX(a);

anfft_a = spectralfftw.partialX(a);
#scipy_a = fft(a, axis=0).real;
#anfft_a = anfft.fft(a.transpose()).real.transpose();


#print scipy_a
#print
#print anfft_a - scipy_a;
