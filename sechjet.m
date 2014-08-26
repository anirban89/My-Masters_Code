clear all
close all
clc

betamin = -3;
betamax = 2;
nbeta = 10;
dbeta = (betamax - betamin)/nbeta;
for j1 = 1:nbeta
    beta(j1) = betamin + (j1-1)*dbeta;
    kmax = 5;
    kmin = 0;
    nk = 10;
    dk = (kmax - kmin)/nk;
    for j = 1:nk
        k(j,j1) = kmin + (j - 1)*dk;
        n = 200;
        yupper = 10;
        ylower = -10;
        d = (yupper - ylower)/n;
% beta = 0.5;
tic
        for i = 1:n+1
            y(i,j,j1) = ylower + d*(i-1);
            U(i,j,j1) = (sech(y(i,j,j1)))^2;
            U11(i,j,j1) = 2*(cosh(2*y(i,j,j1)) - 2)*(sech(y(i,j,j1)))^4;

%     U(i) = (sech(2*y(i) - 1)).^2;
%     U11(i) = -8*(sech(1-2*y(i))).^2 .* ((sech(1-2*y(i))).^2 - 2*(tanh(1-2*y(i))).^2);

            if i > 1 & i < n+1
                A1(i,i,j,j1) = -2 / d^2 - k(j,j1)^2;
                A0(i,i,j,j1) = -(U11(i,j,j1) - beta(j1))*k(j,j1) - U(i,j,j1)* k(j,j1)^3 - 2*U(i,j,j1)* k(j,j1)/ d^2;
                A1(i,i-1,j,j1) = 1/ d^2;
                A0(i,i-1,j,j1) = U(i,j,j1) * k(j,j1) / d^2;
                A1(i,i+1,j,j1) = 1/ d^2;
                A0(i,i+1,j,j1) = U(i,j,j1) * k(j,j1) / d^2;
            else
                A1(i,i,j,j1) = 0;
                A0(i,i,j,j1) = 1;
            end
        end
        [svect(:,:,j,j1),sig(:,:,j,j1)] = eig(A0(:,:,j,j1),A1(:,:,j,j1));
        sig(:,j,j1) = diag(sig(:,:,j,j1));
        sigi(:,j,j1) = imag(sig(:,j,j1));
% plot(y,U,y,svect(:,1),'.r',y,U11,'.k')
    end
    toc
    
% figure, plot(k,sigi,'.r')
end
contour(beta,k,sigi)
