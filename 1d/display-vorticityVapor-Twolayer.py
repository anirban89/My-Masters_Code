from spectral import *;
from integrator import *;
from pylab import *;
import scipy;
import pylab as m;
import time;
import numpy as np;
import sys;
import pylab;


if (len(sys.argv) < 2):
        print 'please give value of kappa as argument';
        exit();

kappa = float(sys.argv[1]);
print kappa;

if (len(sys.argv) == 3):
        M = int(sys.argv[2]);
else:
        M = 10000;

print M;
N = 512;
#initPath = '/media/FreeAgent Drive/Vapor-sameInitialRH/';
initPath = '/home/joy/';
#initPath = '/Data/';
vorticityPath1 = 'TwoLayer'+str(kappa)+'/eta1_data_'
vorticityPath2 = 'TwoLayer'+str(kappa)+'/eta2_data_'
vaporPath = 'TwoLayer'+str(kappa)+'/vapor_data_'
figurePath = 'Vapor2/figure_'

k = arange(0,N/3);
def deriveVelocities(omega1, omega2):
    '''
    omega1 = laplacian(psi1) + (psi2-psi1)
    omega2 = laplacian(psi2) + (psi1-psi2)
    '''
    sumPsi = InvLaplacian(omega1+omega2);
    diffPsi = (omega2 - omega1)/2.0;

    psi1 = (sumPsi + diffPsi)/2.0;
    psi2 = (sumPsi - diffPsi)/2.0;

    return(psi1,psi2);


def calcPower(field):
    [x,y] = shape(field);
    power = zeros(y);

    for i in arange(y):
        val = fft(field[:,i]);
        power = power + abs(val*val.conj());


    power = power/y;
    return log(power[0:N/3]);

def calcVaporPower(field):
    [x,y] = shape(field);
    power = zeros(y);

    for i in arange(y):
        val = fft(field[:,i]);
        power = power + abs(val*val.conj());

    power = power/y;
    return log(power[1:N/2]);


# making my own colormap
cdict = {
    'red'  :  ((0., 1., 1.), (0.5, 0.5, 0.5), (1., 0, 0)),
    'green':  ((0., 1., 1.), (0.5, 0.5, 0.5), (1., 0, 0)),
    'blue' :  ((0., 1., 1.), (0.5, 0.5, 0.5), (1., 0., 0.))
}
my_cmap = m.matplotlib.colors.LinearSegmentedColormap('my_colormap', cdict, 1024)
wavenumber = arange(0,N/3);
wavenumberVapor = arange(0,N/2);
hold(True);
forcing1 = zeros((M,N/3));
forcing2 = zeros((M,N/3));
vaporSpectrum = zeros((M,N/3));
#columnWidth = 237.55974;
columnWidth = 1120/3
inchesToPt = 1.0/72.27;
goldenMean = (sqrt(5)-1.0)/2.0;
figWidth = columnWidth*inchesToPt;
figHeight = figWidth;
figSize = [figWidth, figHeight];

params = {'backend': 'ps',
          'axes.labelsize': 10,
          'text.fontsize': 10,
          'legend.fontsize': 10,
          'xtick.labelsize': 8,
          'ytick.labelsize': 8,
          'figure.figsize': figSize};
params = {'backend': 'ps',
          'axes.labelsize': 18,
          'text.fontsize': 18,
          'legend.fontsize': 18,
          'xtick.labelsize': 10,
          'ytick.labelsize': 10,
          'figure.figsize': figSize};


pylab.rcParams.update(params);


for i in range(M):

    #pvData = np.load(vorticityPath+str(i)+'.npy').transpose();
    value = 100*i;
    pvData1 = np.load(initPath+vorticityPath1+str(value)+'.npy').transpose();
    pvData2 = np.load(initPath+vorticityPath2+str(value)+'.npy').transpose();
    vaporData = np.load(initPath+vaporPath+str(value)+'.npy').transpose();

    (psi1,psi2) = deriveVelocities(pvData1, pvData2);
    print 'reading',vorticityPath1+str(value)+'.npy';
#    pvData[pvData>0]=1;

#    handle.set_array(pvData.transpose());
#    handle.autoscale();

    vorticity1 = partialX(psi1,2)+partialY(psi1,2);


    figure(dpi=800);

    u = -partialY(psi1-psi2);
    v = partialX(psi1-psi2);
    hist(u);
    print 'calculating wind field';
    print 'average U: ', average(abs(u));
    print 'average V: ', average(abs(v));
    u = -partialY(psi2);
    v = partialX(psi2);
    print 'calculating wind field';
    print 'average U: ', average(abs(u));
    print 'average V: ', average(abs(v));
    savefig(initPath+'/LowerLevelVelocity'+str(value)+'.png');
#show();
    print 'plotting vapor field';
    figure(dpi=800);
    pylab.axes([.1,.1,.85,.85]);
    imshow(vaporData.transpose());
    savefig(initPath+'/Vaporfield'+str(value)+'.eps');


#show();

