from spectral import *;
from integrator import *;
from pylab import *;
import scipy;
import pylab as m;
import time;
import sys;


if (len(sys.argv) < 2):
        print 'please give value of kappa as argument';
        exit();

coeff = float(sys.argv[1]);
print coeff;

if (len(sys.argv) == 3):
        time = int(sys.argv[2]);
else:
        time = 10000;

def calcPower(field):
    y = len(field);
    power = zeros(y);

    val = fft(field);
    power = abs(val*val.conj());


    return log(power[0:N/2]);



N = 128;
vorticityPath = 'EqVaportestHiRes'+str(coeff)+'/eta_data_'
#vaporPath = 'EqVaporSuperHiRes'+str(coeff)+'/vapor_data_'
vaporPath = 'EqVaportest'+str(coeff)+'/vapor_data_'
#initPath = '/media/FreeAgent Drive/Vapor'+str(coeff);
initPath = '/home/joy/';
#initPath = '/Data/';
figurePath = 'Figure'+str(coeff)+'/figure_'

#draw();
out = 0;

q0 = 7;
dy = 2*pi/N;
minQ = q0 - coeff*dy*N/2;
offset = 400;

ratioU = zeros(time-400);
anisotropy = zeros(time-400);
ratioNL = zeros(time-400);
timesUpperLimit = 0;
timesLowerLimit = 0;
upperLimitCrossedPrevTime = 0;
lowerLimitCrossedPrevTime = 0;
timeAboveUpperLimit = 0;
timeAboveLowerLimit = 0;

jetDuration = zeros(100);
numJets = 0;

#ratioV = zeros(time);
for i in range(400, time):

    j = 10*i;

    #pvData = np.load(vorticityPath+str(i)+'.npy').transpose();
    pvData = np.load(initPath+vorticityPath+str(j)+'.npy').transpose();
#    vaporData = np.load(initPath+vaporPath+str(j)+'.npy').transpose();
#    print shape(pvData);
    if(i%100 == 0):
        print >> sys.stderr, 'reading',vorticityPath+str(j)+'.npy';
    psi = InvPotentialVorticity(pvData);
    v = partialX(psi);
    u = -partialY(psi);

#    nonlinearTerm = u*partialX(pvData) + v*partialY(pvData);

#    ratioNL[i] = sum(sum(nonlinearTerm*nonlinearTerm))/sum(sum(v*v));

    #ratio[i] = sum(sum(real(u*u + v*v + 10*psi*psi)))/totalPotentialEnergy;
    anisotropy[i-400] = sum(sum(real(u*u))- sum(real(v*v)))/sum(sum(real(u*u + v*v)));
    ratioU[i-400] = sum(sum(real(u*u)))/sum(sum(real(u*u + v*v)));
    if(ratioU[i-400] > 0.8):
        timesUpperLimit = timesUpperLimit+1;
        if(upperLimitCrossedPrevTime):
            timeAboveUpperLimit = timeAboveUpperLimit+1;

        upperLimitCrossedPrevTime = 1;
    else:

        if(upperLimitCrossedPrevTime == 1):
            upperLimitCrossedPrevTime = 0;
            print 'Jet duration: ',str(timeAboveUpperLimit+1);
            jetDuration[numJets] = timeAboveUpperLimit+1;
            numJets = numJets+1;
            timeAboveUpperLimit = 0;


    if(ratioU[i-400] < 0.2):
        timesLowerLimit = timesLowerLimit+1;
#    ratioV[i] = sum(sum(real(v*v)))/sum(sum(real(u*u + v*v)));
    #ratio[i] = sum(sum(real(10*psi*psi)))/totalPotentialEnergy;

#    ratio[i] = sum(sum(real(u*u + v*v)));

#    out = raw_input();
#    out = int(out);

    #scipy.misc.toimage(pvData, cmin=0, cmax=255).save(figurePath+str(i)+".png");
    #imsave(figurePath+str(i)+".png",pvData);

plot(ratioU, label='u^2/KE');
title('ratio of zonal KE to total KE for kappa = '+str(coeff));
hold(True);
#plot(ratioV, label = 'v^2/KE');
legend();
savefig(initPath+'/ZonalVel'+str(int(coeff*100))+'.png');
print 'Upper Limit Crossed '+str(timesUpperLimit)+' times';
print 'Lower Limit Crossed '+str(timesLowerLimit)+' times';

figure();
print >> sys.stderr, numJets;

if(numJets > 0):
    bins = [0,5,10,15,20,25,30,35,40,45,50,55,60,65];
    (hist, bins) = histogram(jetDuration[0:numJets],bins);
    bins = (bins[0:len(bins)-1] + bins[1:len(bins)])/2.0;
    print (bins), len(hist);
    bar(bins*10, hist);
    yticks(range(0,50,5));
    title('\nMax Jet Width ='+str(amax(jetDuration[0:numJets])*10)+' time steps'+'\n Number of Jets: '+str(numJets));
    savefig(initPath+'/JetDuration'+str(int(coeff*100))+'.png');
    figure();
#psd(ratioU - average(ratioU), time-400, 0.1);
plot(calcPower(ratioU - average(ratioU)));
savefig(initPath+'/PSD'+str(coeff)+'.png');
#show();
