function [dnew,y]=jd1(yb,d);

[yb,d]=chebd(200,0,1);
% y0=30;sp=10;                            % High Sp to ensure stability of d4  -- Verify

y=atanh(yb);

% y=sp*(1+yb)./(1-yb+2*sp/y0);			%Map y [0 y0] to [-1 1]		--- NASA 1992 Mapping
% %plot(yb,y,'o');
% %hold all;
% j=(1-yb+2*sp/y0).^2/(2*sp*(1+sp/y0));		% Jacobian for the transformation	
%plot(yb,j,'b');

j = 1 - yb.^2;
j = diag(j);
dnew = j*d;

