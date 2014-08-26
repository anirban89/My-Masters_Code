from pylab import *;

N = 24;

beta = 1;
y = linspace(-5,5,11);
k = linspace(0,10,11);

matrix = zeros((N,N),dtype=complex);
'''
matrix[0,1] = beta*y[0];
matrix[0,2] = 1j*k[1];

matrix[1,0] = beta*y[0];
matrix[1,5] = 0.5;

matrix[2,0] = 1j*k[1];
matrix[2,4] = 0.5;

matrix[N-3,N-2] = beta*y[(N-3)/3];
matrix[N-3,N-1] = 1j*k[1];

matrix[N-2,N-3] = beta*y[(N-3)/3];
matrix[N-2,N-1] = -0.5;

matrix[N-1,N-3] = 1j*k[1];
matrix[N-1,N-2] = -0.5;
'''
#for row in arange(1,(N-3)/3):
for row in arange(0,(N-3)/3):
    matrix[3*row, 3*row+1] = beta*y[row];
    matrix[3*row, 3*row+2] = 1j*k[1];

    matrix[3*row+1,3*row] = beta*y[row];
    if(row>0):
        matrix[3*row+1,3*row-1] = -0.5;
    matrix[3*row+1,3*row+5] = 0.5;

    matrix[3*row+2,3*row] = 1j*k[1];
    if(row>0):
        matrix[3*row+2,3*row-2] = -0.5;
    matrix[3*row+2,3*row+4] = 0.5;


#print real(matrix)+imag(matrix);

(V,M) = eig(matrix);
#print V
plot(real(V));
hold(True);
plot(imag(V));
show();
