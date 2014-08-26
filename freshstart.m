clear all
close all
clc

% for j = 1:10
% k(j) = j - 1;
% d = .2;
% beta = 11;
% for i = 1:4
%     y(i,j) = 0.2 * i;
%     A1(i,i,j) = -2 / d^2 - k(j)^2;
%     A0(i,i,j) = -k(j) * (pi)^2 * sin(pi * y(i,j)) + sin(pi * y(i,j)) * k(j) * (2/ d^2 + k(j)^2) - k(j) * beta; 
%     if i > 1
%         A1(i,i-1,j) = 1/ d^2;
%         A0(i,i-1,j) = - sin(pi * y(i,j)) * k(j) * (1/ d^2);
%     end
%     if i < 4
%         A1(i,i+1,j) = 1/ d^2;
%         A0(i,i+1,j) = - sin(pi * y(i,j)) * k(j) * (1/ d^2);
%     end
%     
% end
% [svect(i), sig(i)] = polyeig (A0(i),A1(i));
% sigr = real(sig);
% end


kmax = 6;
kmin = 0;
nk = 100;
dk = (kmax - kmin)/nk;
n = 200;
yupper = 10;
ylower = -10;
d = (yupper - ylower)/n;
beta = -1.99;
tic
for j = 1:nk
    k(j) = kmin + (j - 1)*dk;

for i = 1:n+1
     y(i,j) = ylower + d*(i-1);
%      U(i) = sin(pi*y(i));
%      U11(i) = - pi^2 * sin(pi*y(i));
%  
          U(i,j) = (sech(y(i,j)))^2;
          U11(i,j) = 2*(cosh(2*y(i,j)) - 2)*(sech(y(i,j)))^4;

%     U(i) = (sech(2*y(i) - 1)).^2;
%     U11(i) = -8*(sech(1-2*y(i))).^2 .* ((sech(1-2*y(i))).^2 - 2*(tanh(1-2*y(i))).^2);

    if i > 1 & i < n+1
         A1(i,i,j) = -2 / d^2 - k(j)^2;
         A0(i,i,j) = -(U11(i,j) - beta)*k(j) - U(i,j)* k(j)^3 - 2*U(i,j)* k(j)/ d^2;
         A1(i,i-1,j) = 1/ d^2;
         A0(i,i-1,j) = U(i,j) * k(j) / d^2;
         A1(i,i+1,j) = 1/ d^2;
         A0(i,i+1,j) = U(i,j) * k(j) / d^2;
     else
         A1(i,i,j) = 0;
         A0(i,i,j) = 1;
     end
end
[svect(:,:,j),sig(:,:,j)] = eig(A0(:,:,j),A1(:,:,j));
sigma(:,j) = diag(sig(:,:,j));
sigi(:,j) = imag(sigma(:,j));
sigr(:,j) = real(sigma(:,j));
% plot(y,U,y,svect(:,1),'.r',y,U11,'.k')
end
toc


figure, plot(k,sigi,'.r')