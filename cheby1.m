clear all
close all
clc
num = 100;
[r,D,D2,D3,D4] = trafoJet(num,10,20);
W0 = (sech (r - 10)).^2;
% W0 = sin(2*pi*r);
[eval, evec] = rotation(r, D, D2, W0, 1.2);
plot(r,W0,r,D2*W0,'.k',r,evec(:,1),'.r')
