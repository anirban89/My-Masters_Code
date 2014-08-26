from scipy import linalg,pi,array,zeros,arange;
from threadMultiply import *;


class coupledDerivative:

    def __init__(self, xLen, yLen, xSize=None, ySize=None):
        """
            xLen,yLen = number of grid points required
            xSize, ySize = physical size of domain

            The init calculates the required matrices and stores them
            for future calculations.


        """
      
        if(xSize == None):
            xSize = 6*pi;

        if(ySize == None):
            ySize = 2*pi;

        self.xSize = xSize;
        self.xLen = xLen;
        self.yBufferBC = 20;
        self.ySize = ySize;
        self.yLen = yLen + self.yBufferBC;#Need to add boundary conditions. keep
                             #place for that
        self.queue = Queue.Queue();
        self.numThreads = 4;
         
        h = xSize/xLen;
        self.xh = xSize/xLen;
        matrix = zeros((2*xLen,2*xLen));
        row1 = array([51, 9*h, 108, 0, 51, -9*h]);
        row2 = array([-138, -18*h, 0, 108*h, 138, -18*h]);
        index = arange(0,6);

        for i in arange(0,2*xLen,2):
            newIndex = (index-2+i)%(2*xLen);
            matrix[i,newIndex] = row1[:];
            matrix[i+1,newIndex] = row2[:];

        self.xMatrix =  linalg.inv(matrix);
        print 'Stored X Matrix';

        yLen = yLen + self.yBufferBC;#We need to explicitly add BCs

        h = 1;
        h = ySize/(yLen - self.yBufferBC);
        self.hy = ySize/(yLen - self.yBufferBC);
        matrix = zeros((2*yLen,2*yLen));
        row1 = array([51, 9*h, 108, 0, 51, -9*h]);
        row2 = array([-138, -18*h, 0, 108*h, 138, -18*h]);
        index = arange(0,6);

        for i in arange(0,2*yLen,2):
            newIndex = (index+i)%(2*yLen);
            matrix[i,newIndex] = row1[:];
            matrix[i+1,newIndex] = row2[:];

        #print matrix;
        self.yMatrix =  linalg.inv(matrix);
        print 'Stored Y Matrix';

        #Spawn numThreads threads to do the calculations.

        for i in arange(self.numThreads):
            t = threadMultiply(self.queue, self.xLen,
                    self.xMatrix, self.yMatrix);
            t.setDaemon(True);
            t.start();

        print 'Spawned '+str(self.numThreads)+' threads';

    def partialX(self, inp, returnDoublePrime=False):

        i = arange(self.xLen);
        inMinusOne = inp[i-1,:];
        inMinusTwo = inp[i-2,:];
        inPlusOne = inp[(i+1)%(self.xLen),:];
        inPlusTwo = inp[(i+2)%(self.xLen),:];


        indexEven = arange(0,2*self.xLen,2);
        indexOdd = arange(1,2*self.xLen,2);
        inputData = zeros((2*self.xLen,self.yLen-self.yBufferBC));
        outData = zeros((2*self.xLen,self.yLen-self.yBufferBC));

        #Get data in required form for calculations

        inputData[indexEven,:] = 107*(inPlusOne-inMinusOne) - (inPlusTwo-inMinusTwo);
        inputData[indexOdd,:] = 352*(inPlusOne+inMinusOne) - (inPlusTwo+inMinusTwo) - 702*inp;
        inputData = inputData/(self.xSize/self.xLen);


        for i in arange(self.numThreads):
            self.queue.put((arange((self.yLen-self.yBufferBC)/self.numThreads)+
                i*(self.yLen-self.yBufferBC)/self.numThreads, inputData, outData, 'X'));

        #Wait for calculations to finish
        self.queue.join();

        if(returnDoublePrime == True):
            return (outData[indexEven,:], outData[indexOdd,:]);
        else:
            return (outData[indexEven,:]);

    def partialY(self, inp, bc, returnDoublePrime=False):

        i = arange(self.yLen);
        inp = bc.applyBoundaryNS(inp);

        inMinusOne = inp[:,i-1];
        inMinusTwo = inp[:,i-2];
        inPlusOne = inp[:,(i+1)%(self.yLen)];
        inPlusTwo = inp[:,(i+2)%(self.yLen)];


        indexEven = arange(0,2*self.yLen,2);
        indexOdd = arange(1,2*self.yLen,2);
        inputData = zeros((self.xLen,2*self.yLen));
        outData = zeros((self.xLen,2*self.yLen));


        #Get data in required form for calculations
        
        inputData[:,indexEven] = 107*(inPlusOne-inMinusOne) - (inPlusTwo-inMinusTwo);
        inputData[:,indexOdd] = 352*(inPlusOne+inMinusOne) - (inPlusTwo+inMinusTwo) - 702*inp;
        inputData = inputData/(self.ySize/(self.yLen-self.yBufferBC));


        for i in arange(self.numThreads):
            self.queue.put((arange(self.xLen/self.numThreads)+i*self.xLen/self.numThreads,
            inputData.transpose(), outData.transpose(), 'Y'));

        #Wait for calculations to finish
        self.queue.join();


        if(returnDoublePrime == True):
            return \
            (outData[:,indexEven[(self.yBufferBC/2+1):(len(indexEven)+1-self.yBufferBC/2)]],
            outData[:,indexOdd[(self.yBufferBC/2+1):(len(indexOdd)+1-self.yBufferBC/2)]]);
        else:
            return (outData[:,indexEven[self.yBufferBC/2:len(indexEven)-self.yBufferBC/2]]);



def main():
    
    a = coupledDerivative(8,2);

#main();
