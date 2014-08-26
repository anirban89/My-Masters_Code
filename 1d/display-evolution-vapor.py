from spectral import *;
from integrator import *;
from pylab import *;
import scipy;
import pylab as m;
import time;
import sys;

N = 128;

if(len(sys.argv) < 2):
    exit;

coeff = float(sys.argv[1]);
print coeff;

if(len(sys.argv) == 3):
    time = int(sys.argv[2]);
else:
    time = 10000;
#initPath = '/media/FreeAgent Drive/Vapor-sameInitialRH/';
vorticityPath = 'EqVaportestHiRes'+str(coeff)+'/eta_data_'
initPath = '/home/joy/';
#initPath = '/Data/';
#vaporPath = 'EqVaporSuperHiRes'+str(coeff)+'/vapor_data_'
#vaporPath = 'EqVaportestHiResBoundaryForcing'+str(coeff)+'/vapor_data_'
#vaporPath = 'EqVaportestLoRes'+str(coeff)+'/vapor_data_'
vaporPath = 'EqVaportestHiRes'+str(coeff)+'/vapor_data_'
deltaPath = 'Vapor'+str(coeff)+'/delta_data_'
figurePath = 'Figure'+str(coeff)+'/figure_'

q0 = 7;
dy = 2*pi/N;

saturatedFieldG = zeros((N,N));

for i in arange(N/2):
    saturatedFieldG[:,i] = q0 - coeff*dy*(N/2-i);
    saturatedFieldG[:,N-i-1] = q0 - coeff*dy*(N/2-i);

ion()
# making my own colormap
cdict = {
    'red'  :  ((0., 1., 1.), (0.5, 0.5, 0.5), (1., 0, 0)),
    'green':  ((0., 1., 1.), (0.5, 0.5, 0.5), (1., 0, 0)),
    'blue' :  ((0., 1., 1.), (0.5, 0.5, 0.5), (1., 0., 0.))
}
my_cmap = m.matplotlib.colors.LinearSegmentedColormap('my_colormap', cdict, 1024)

handle = imshow(zeros((N,N)));
colorbar();
draw();
out = 0;

offset = 000;

vaporSeries = zeros(time-offset);
subSat = zeros((time-offset));
superSat = zeros((time-offset));
justSat = zeros((time-offset));


for i in range(offset,time):

    j = 10*i;

#    pvData = np.load(initPath+vorticityPath+str(j)+'.npy').transpose();
    vData = np.load(initPath+vaporPath+str(j)+'.npy').transpose();
#    deltaData = np.load(initPath+deltaPath+str(i)+'.npy').transpose();
    print 'reading',vaporPath+str(j)+'.npy';
    print shape(vData);
    vData = vData;
    saturatedField = saturatedFieldG;
    print shape(saturatedField);
    vaporField = vData+saturatedField;
    superSat[i-offset] =  sum(vaporField[vData>0]);
    subSat[i-offset] =  -sum(vaporField[vData<0]);
    justSat[i-offset] = sum(vaporField[vData==0]);
#    rhData = (vaporField)/saturatedField;
#    vaporSeries[i] = sum(sum(vData*vData));


#    psi = InvPotentialVorticity(pvData);

#    u = partialX(psi)

    handle.set_array(vData.transpose());
    handle.autoscale();
    draw();
   # time.sleep(0.05);
#    out = raw_input();
#    out = int(out);

    #scipy.misc.toimage(pvData, cmin=0, cmax=255).save(figurePath+str(i)+".png");
    #imsave(figurePath+str(i)+".png",pvData);

exit();
timeAxis = arange(offset,time);
print len(timeAxis);
print len(subSat);
hold(True);
figure(figsize=(16,10));
sumSat = justSat+subSat+superSat;
plot(timeAxis, subSat/amax(subSat), label='SubSat');
plot(timeAxis, superSat/amax(superSat), label='SuperSat');
#plot(timeAxis, justSat, label='JustSat');
legend(loc='lower left');
#plot(timeAxis, vaporSeries/max(vaporSeries));
savefig(initPath+'VaporEvolution'+str(int(100*coeff))+".png");

figure(figsize=(16,10));
spectrum = fft(superSat);
plot(log(spectrum*conj(spectrum))[0:len(spectrum)/2]);
savefig(initPath+'VaporTimeSeriesSpectrum'+str(int(100*coeff))+".png");
#show();


