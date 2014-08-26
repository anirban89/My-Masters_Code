from pylab import *;
from finiteDifference import *;
import time;

N = 1024;
h = 2*pi/N;
#h = 1;

i = arange(N);

inp = sin(36*pi*i/N);
inpPrime = 18*cos(36*pi*i/N);
inpDoublePrime = -324*sin(36*pi*i/N);

inMinusOne = inp[i-1];
inMinusTwo = inp[i-2];
inPlusOne = inp[(i+1)%N];
inPlusTwo = inp[(i+2)%N];


indexEven = arange(0,2*N,2);
indexOdd = arange(1,2*N,2);
inputData = zeros(2*N);

inputData[indexEven] = 107*(inPlusOne-inMinusOne) - (inPlusTwo-inMinusTwo);
inputData[indexOdd] = 352*(inPlusOne+inMinusOne) - (inPlusTwo+inMinusTwo) - 702*inp;
inputData = inputData/h;

matrix = zeros((2*N,2*N));
row1 = array([51, 9*h, 108, 0, 51, -9*h]);
row2 = array([-138, -18*h, 0, 108*h, 138, -18*h]);
index = arange(0,6);

for i in arange(0,2*N,2):
    newIndex = (index-2+i)%(2*N);
    matrix[i,newIndex] = row1[:];
    matrix[i+1,newIndex] = row2[:];

print matrix;
inverse =  inv(matrix);
print shape(inputData);
print shape(inverse);

a = zeros((N*2,N));
for i in arange(N):
    a[:,i] = inputData;

start = time.time();
out =  dot(inverse,a);

print 'Time: ', time.time()-start;
print(shape(out));
print(shape(inpPrime));

plot(inpPrime - out[indexEven,0]);
figure();
plot(inpDoublePrime - out[indexOdd,0]);

figure();
newInp = zeros((N));
newInp[:] = inp[:];
plot(inpPrime - boydPartialX(newInp));
figure();
plot(inpDoublePrime - boydPartialX(boydPartialX(newInp)));
show();


