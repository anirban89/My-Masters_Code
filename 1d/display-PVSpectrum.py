from spectral import *;
from integrator import *;
from pylab import *;
import scipy;
import pylab as m;
import time;
import numpy as np;

N = 530;
kappa = 0.1;
initPath = '/media/FreeAgent Drive/Vapor'+str(coeff);
vorticityPath = '/eta_data_'
vaporPath = 'Vapor'+str(kappa)+'/vapor_data_'
figurePath = 'Vapor2/figure_'

def calcPower(field):
    [x,y] = shape(field);
    power = zeros(y);

    for i in arange(x):
        power = power + abs(fft(field[i,:]));

    power = power/x;
    return log(power);

ion()
# making my own colormap
cdict = {
    'red'  :  ((0., 1., 1.), (0.5, 0.5, 0.5), (1., 0, 0)),
    'green':  ((0., 1., 1.), (0.5, 0.5, 0.5), (1., 0, 0)),
    'blue' :  ((0., 1., 1.), (0.5, 0.5, 0.5), (1., 0., 0.))
}
my_cmap = m.matplotlib.colors.LinearSegmentedColormap('my_colormap', cdict, 1024)
forcing = zeros((N))
corr = zeros((N));

handle, = plot(zeros(256));
draw();

for i in range(0,510,100):

    #pvData = np.load(vorticityPath+str(i)+'.npy').transpose();
    pvData = np.load(initPath+vorticityPath+str(i)+'.npy').transpose();
    print shape(pvData);
    print 'reading',vorticityPath+str(i)+'.npy';
    spec = calcPower(pvData);

    figure();
    plot(spec);
    draw();
    out = raw_input();
    time.sleep(0.5);

