from pylab import *

from boundary import *


bc = boundary();

a = rand(4,4);
print a
print '-------------------------------------------------'
print
print 'applying north south conditions'

print bc.applyBoundaryNS(a);

print '-------------------------------------------------'
print
print 'applying east west conditions'
print bc.applyBoundaryEW(a);

print '-------------------------------------------------'
print
print 'applying complete boundary'
print bc.applyCompleteBoundary(a);

a = rand(2,5);

print
print 'New a:'
print a.transpose();
bc2 = boundary(yBcType='Symmetric');
print '-------------------------------------------------'
print
print 'applying Symmetric boundary'
print bc2.applyBoundaryNS(a).transpose();
bc3 = boundary(yBcType='AntiSymmetric');
print '-------------------------------------------------'
print
print 'applying AntiSymmetric boundary'
print bc3.applyBoundaryNS(a).transpose();
