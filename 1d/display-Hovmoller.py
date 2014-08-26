from spectral import *;
from integrator import *;
from pylab import *;
import scipy;
import pylab as m;
import time;
import numpy as np;
import sys;
import pylab;
from os.path import expanduser;

home = expanduser("~");

if(len(sys.argv) > 1):
    if(len(sys.argv) == 4):
    
        alphaVapor = float(sys.argv[1]);
        betaVapor = float(sys.argv[2]);
        rangeVals = int(sys.argv[3]);
        print "alphaQ = ", str(alphaVapor);
        print "betaQ = ", str(betaVapor);
else:
    print 'Insufficient arguments';
    exit();

N = 256;
h = 2*pi/N;

tn = 0;
sec = 1;

initPath = home+'/LinearRuns/TwoLayerLinear';

if alphaVapor<0:
    initPath = home+'/LinearRuns/TwoLayerLinearNegative';

path = initPath+str(int(10*abs(alphaVapor)))+str(int(10*betaVapor));
figPath = home+'/lowResRuns/'+str(int(10*abs(alphaVapor)))+str(int(10*betaVapor));


dx = 2*pi/(N-1);
dy = 2*pi/(N-1);
qs = zeros((N,N),dtype=double);
q0 = 7;

betaMax = q0/(2*pi);
alphaMax = q0/(2*pi);

beta = betaVapor*betaMax;
alpha = (alphaVapor*alphaMax);

if(alphaVapor > 0):
    for i in arange(N):
        for j in arange(N):
            qs[N-i-1,j] = q0 - alpha*dx*(i) - beta*dy*j;
elif(alphaVapor < 0):
    for i in arange(N):
        for j in arange(N):
            qs[i,j] = q0 - alpha*dx*(i) - beta*dy*j;
else:
    for i in arange(N):
        for j in arange(N):
            qs[i,j] = q0 - beta*dy*j;



#initPath = '/media/FreeAgent Drive/Vapor-sameInitialRH/';
#initPath = '/Data/';
vorticityPath1 = '/eta1_data_'
vorticityPath2 = '/eta2_data_'
vaporPath = '/vapor_data_'

k = arange(0,N/3);
def deriveVelocities(omega1, omega2):
    '''
    omega1 = laplacian(psi1) + (psi2-psi1)
    omega2 = laplacian(psi2) + (psi1-psi2)
    '''
    sumPsi = InvLaplacian(omega1+omega2);
    diffPsi = InvPotentialVorticityTwoLayer(omega2 - omega1);

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

vaporHovmoller = zeros((rangeVals,N));
pvHovmoller = zeros((rangeVals,N));
condensation = zeros(rangeVals);
vaporVals = zeros(rangeVals);
figure();
for i in arange(rangeVals):

    j = (i+1)*10;

    try:
    #pvData = np.load(vorticityPath+str(i)+'.npy').transpose();
        pvData1 = np.load(path+vorticityPath1+str(j)+'.npy').transpose();
        pvData2 = np.load(path+vorticityPath2+str(j)+'.npy').transpose();
        vaporData = np.load(path+vaporPath+str(j)+'.npy');
        vaporData[vaporData<0] = 0;
        print 'reading time step:', str(j);
        print shape((average(vaporData,axis=0)));
        vaporHovmoller[i,:] = average(vaporData,axis=0);
        deviation = max(vaporHovmoller[i,:]) - min(vaporHovmoller[i,:]);
        print average(vaporHovmoller[i,:]);
        normalised = (vaporHovmoller[i,:]/abs(deviation));
        normalised = normalised-average(normalised);
        vaporHovmoller[i,:] = normalised;

#        psi1,psi2 = deriveVelocities(pvData1, pvData2);
#        baroclinic = psi2;
#        pvHovmoller[i,:] = average(baroclinic,axis=1);
#        pvHovmoller[i,:] = average(pvData1,axis=1);
#        deviation = max(pvHovmoller[i,:]) - min(pvHovmoller[i,:]);
#        print 'U Vel: ', amax(fabs(partialY(psi2)));
#        print 'V Vel: ', amax(fabs(partialX(psi2)));
#        normalised = (pvHovmoller[i,:]/abs(deviation));
#        normalised = normalised-average(normalised);
#        pvHovmoller[i,:] = normalised;

        heaviside = zeros((N,N));
        heaviside[vaporData>0] = 1;
        condensation[i] = sum(heaviside);
#vaporVals[i] = sum(vaporData[vaporData>0]);
        vaporVals[i] = sum(qs+vaporData);
    except:
        vaporVals[i] = vaporVals[i-1];
        pvHovmoller[i,:] = pvHovmoller[i-1,:];
        pass;

    '''plot(normalised);
    draw();

    a = raw_input();
'''
#    (psi1,psi2) = deriveVelocities(pvData1, pvData2);
#    print 'reading',vorticityPath1+str(value)+'.npy';
#    pvData[pvData>0]=1;

#    handle.set_array(pvData.transpose());
#    handle.autoscale();

#    vorticity1 = partialX(psi1,2)+partialY(psi1,2);

#    out1 = calcPower(vorticity1);
#    out2 = calcPower(vorticity1.transpose());
#    out3 = calcPower(vaporData);
    #draw();
#    forcing1[i,:] =  out1;
#    forcing2[i,:] =  out2;
#    vaporSpectrum[i,:] =  out3*k;
    #out = raw_input();
    #exit()
#    out = int(out);

    #scipy.misc.toimage(pvData, cmin=0, cmax=255).save(figurePath+str(i)+".png");
    #imsave(figurePath+str(i)+".png",pvData);
'''
columnWidth = 237.55974;
inchesToPt = 1.0/72.27;
goldenMean = (sqrt(5)-1.0)/2.0;
figWidth = columnWidth*inchesToPt;
figHeight = figWidth*goldenMean;
figSize = [figWidth, figHeight];

params = {'backend': 'ps',
          'axes.labelsize': 8,
          'text.fontsize': 8,
          'legend.fontsize': 10,
          'xtick.labelsize': 8,
          'ytick.labelsize': 8,
          'figure.figsize': figSize};

pylab.rcParams.update(params);
'''
#pylab.axes([.1,.2,.85,.75]);
subplot(211);
plot((vaporVals.transpose()))

subplot(212);
imshow(vaporHovmoller.transpose(),aspect='auto',extent=(0,rangeVals,0,2*pi));
#colorbar();
'''
vaporVals[1:-2] = (vaporVals[2:-1]-vaporVals[0:-3])/2.0;
vaporVals[0] = vaporVals[-1] = 0;
subplot(313);
imshow(pvHovmoller.transpose(),aspect='auto');
'''

savefig(figPath+'VaporHovmoller.jpg',dpi=400);
show();

