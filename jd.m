function [dnew,y]=jd(yb,d);

[yb,d]=cheb1(161);
% y0=20;sp=10;                            % High Sp to ensure stability of d4  -- Verify
% L = 30; 
% L1 = -L ; L2 = -1/L;
% y=atanh(yb);
% y = (L1+yb)./(L2+yb);
% y = y0 * yb;

n1 = 161;
p = 0.015;
q = 0.12;
a1 = p + q/2;
b1 = 8;
y0 = (0.5/b1)*log( (1+(exp(b1) - 1)*a1)/(1+(exp(-b1) - 1)*a1));
Y11 = sinh(b1*y0) ;
yl = 10;
ihalf = (n1 - 1)/2 +1;
y = zeros(n1,1);
for i = 1:ihalf
    ymid = yb(i);
y(i) = yl*(a1/Y11).*(sinh((ymid - y0)*b1) + Y11);
y12 = (ymid - y0)*b1;
j(i,i) = Y11/(a1*b1*yl.*(cosh(y12)));
end
for i = ihalf+1:n1
    isym = n1+1-i;
    y(i) = -y(isym);
    j(i,i) = j(isym,isym);
end



% y=sp*(1+yb)./(1-yb+2*sp/y0);			%Map y [0 y0] to [-1 1]		--- NASA 1992 Mapping
% %plot(yb,y,'o');
% %hold all;
% j=(1-yb+2*sp/y0).^2/(2*sp*(1+sp/y0));		% Jacobian for the transformation	

% j = 1/y0;
% j = (L2 + yb).^2 /(L2 - L1);
% j = 1 - yb.^2;
% j = diag(j);
dnew = j*d;

