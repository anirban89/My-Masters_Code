function [r,D] = chebd(n,rc,rmax)

% define transformation xi->y
xi = cos( (0:n-1)*pi/(n-1) );
%r = rc*(1-xi)./(1-xi.^2+2*rc/rmax);
r = rc*(1-xi)./(1+xi+2*rc/rmax);
%r = rc*(1+xi)./(1-xi+2*rc/rmax);
r = reshape(r,length(r),1);  % make sure r is a column vector

% build derivative matrices D and DD
% such that F'(i) = D(i,:)*F(:) and F''(i) = DD(i,:)*F(:)
% with F(i) = F(r(i)) in r(i) as defined above.
A=zeros(n);
D=zeros(n);
C = ones(1,n);  C(1) = 2;  C(n) = 2;
for i=1:n
    for j=1:n
        if (i~=j)
            A(i,j) = C(i)*(-1)^(i+j) / (C(j)*(xi(i)-xi(j)));
        elseif (i==1)
            A(i,i) = (2*(n-1)^2 + 1)/6;
        elseif (i==n)
            A(i,i) = -(2*(n-1)^2 + 1)/6;
        else
            A(i,i) = -xi(i) / (2*(1-xi(i)^2));
        end
    end
end
% S = dxi/dr at r(i)
S = zeros(size(r));
%S(1) = -(2*(1 + rc/rmax))/rc;
for i=1:n
%S(i) = (rc^2 - 2*rc*r(i))/(2*r(i)^3*sqrt( ((rc-2*r(i))/r(i))^2 + 8*rc/rmax )) -rc/(2*r(i)^2);
S(i) = -(2*rc*(1 + rc/rmax))/((r(i) + rc)^2);
%S(i) = (2*(1+rc/rmax)/rc)*((r(i)/(1+rc))^2);
end

for i=1:n
    for j=1:n
        D(i,j) = S(i)*A(i,j);
    end
end

% D2 = D*D;
% D3 = D*D2;
% D4 = D*D3;
