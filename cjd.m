function [dnew,y]=sjd(sp,N)



clear all
N=150;
[yb,d]=chebdif(N,1);

y=atanh(yb);
%y2=asinh(30*yb);
y3=30*yb/(1-yb.^2).^0.5;

sp=10;
y0=30;					% High Sp to ensure stability of d4  -- Verify
y4=sp*(1+yb)./(1-yb+2*sp/y0);			%Map y [0 y0] to [-1 1]		--- NASA 1992 Mapping
%plot(yb,y,'o');
%hold all;
j=(1-yb+2*sp/y0).^2/(2*sp*(1+sp/y0));		% Jacobian for the transformation	
%plot(yb,j,'b');
j = diag(j);
dnew = j*d;

