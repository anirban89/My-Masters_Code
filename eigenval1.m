clear all
close all
clc



kmax = 3;
kmin = 0;
nk = 60;
dk = (kmax - kmin)/nk;
k = kmin:dk:kmax;
n = 200;
yupper = 10;
ylower = -10;
d = (yupper - ylower)/n;
y = ylower:d:yupper;
U = (sech(y)).^2;
U11 = 2*(cosh(2*y) - 2).*(sech(y)).^4;
blower = -2.5;
bupper = 1;
nbeta = 70;
dbeta = (bupper - blower)/nbeta ;
beta = blower:dbeta:bupper;
% tic
for m = 1:nbeta +1
    
for j = 1:nk+1
   

for i = 1:n+1
     

    if i > 1 && i < n+1
         
         A1(i,i) = -2 / d^2 - k(j)^2;
         A0(i,i) = -(U11(i) - beta(m))*k(j) - U(i)* k(j)^3 - 2*U(i)* k(j)/ d^2;
         A1(i,i-1) = 1/ d^2;
         A0(i,i-1) = U(i) * k(j) / d^2;
         A1(i,i+1) = 1/ d^2;
         A0(i,i+1) = U(i) * k(j) / d^2;
    else
         A1(i,i) = 0;
         A0(i,i) = 1;
    end
     
end

[svect,sig] = eig(A0,A1);

sigi = imag(sig);
gr(j,m) = max(max(sigi));
%ind = find(sigi(:,j) > 0 );

clear A0 A1

end
end
contour(beta,k,gr);
xlabel('\beta')
ylabel('k')
title('eigenvalues for U = sech^2(y) finite difference')
colorbar

