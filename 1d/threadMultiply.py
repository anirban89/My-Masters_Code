import threading;
from pylab import *;
import Queue;
import time;

N = 1024;
h = 2*pi/N;
a = zeros((N*2,N));
b = zeros((N*2,N));

class threadMultiply(threading.Thread):

    
    def __init__(self,queue,N,xMatrix,yMatrix):

        self.const = 2;
        threading.Thread.__init__(self);
        self.queue = queue;

        self.xMatrix =  xMatrix;
        self.yMatrix =  yMatrix;
        print 'Stored required matrices';

    def run(self):
    
        while True:

            (index, inputData, outputData, partialAlong) = self.queue.get();
            if(partialAlong == 'X'):\
                matrix = self.xMatrix;
            else:
                matrix = self.yMatrix;

            outputData[:,index] = dot(matrix,
                    inputData[:,index]);

            self.queue.task_done();


def main():

    queue = Queue.Queue();


    matrix = zeros((2*N,2*N));
    row1 = array([51, 9*h, 108, 0, 51, -9*h]);
    row2 = array([-138, -18*h, 0, 108*h, 138, -18*h]);
    index = arange(0,6);

    for i in arange(0,2*N,2):
        newIndex = (index-2+i)%(2*N);
        matrix[i,newIndex] = row1[:];
        matrix[i+1,newIndex] = row2[:];

    matrix = inv(matrix);
    print 'Done inverting Matrix';
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

    for i in arange(N):
        a[:,i] = inputData;


    for i in arange(5):

        t = threadMultiply(queue,a,b,N,matrix);
        t.setDaemon(True);
        t.start();


    start = time.time();
    numThreads = 4;
    for i in arange(numThreads):

        queue.put(arange(N/numThreads)+i*N/numThreads);
        #queue.put(i);

    
    queue.join();



    print 'Time taken: ', time.time()-start;

    plot(b[indexOdd,30]-inpDoublePrime);
    hold(True);
    plot(b[indexEven,30]-inpPrime);
    show();


#main();


            


