from spectral import *;
from integrator import *;
from pylab import *;
import scipy;
import pylab as m;
import time;
import numpy as np;
import sys;
import matplotlib.pyplot as plt;
import pylab;


if (len(sys.argv) < 2):
        print 'please give value of kappa as argument';
        exit();

kappa = float(sys.argv[1]);
print kappa;

if (len(sys.argv) == 3):
        value = int(sys.argv[2]);
else:
        value = 10000;

print value;
N = 512;
#initPath = '/media/FreeAgent Drive/Vapor-sameInitialRH/';
initPath = '/home/joy/';
#initPath = '/Data/';
vorticityPath1 = 'TwoLayer'+str(kappa)+'/eta1_data_'
vorticityPath2 = 'TwoLayer'+str(kappa)+'/eta2_data_'
vaporPath = 'TwoLayer'+str(kappa)+'/vapor_data_'
figurePath = 'Vapor2/figure_'

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

x = np.arange(512)
y = np.arange(512)
X, Y = np.meshgrid(x, y)

columnWidth = 237.55974;
inchesToPt = 1.0/72.27;
goldenMean = (sqrt(5)-1.0)/2.0;
figWidth = columnWidth*inchesToPt;
figHeight = figWidth;
figSize = [figWidth, figHeight];

params = {'backend': 'ps',
          'axes.labelsize': 8,
          'text.fontsize': 8,
          'legend.fontsize': 10,
          'xtick.labelsize': 8,
          'ytick.labelsize': 8,
          'figure.figsize': figSize};

pylab.rcParams.update(params);

figure(dpi=800);
pylab.axes([.1,.1,.85,.85]);


pvData1 = np.load(initPath+vorticityPath1+str(value)+'.npy').transpose();
pvData2 = np.load(initPath+vorticityPath2+str(value)+'.npy').transpose();
vaporData = np.load(initPath+vaporPath+str(value)+'.npy').transpose();

(psi1,psi2) = deriveVelocities(pvData1, pvData2);
print 'reading',vorticityPath1+str(value)+'.npy';


manual_locs = concatenate((linspace(.3,2.4,30),linspace(-2.4,-0.3,30))); 
print manual_locs;
CS = plt.contour(X,Y,psi1,levels=manual_locs);
plt.clabel(CS);

savefig(initPath+'/HighVaporContours.eps');
plt.show();
