function [x, D] = cheb1(N)

x = cos( (0:N-1)*pi/(N-1) );

x = reshape(x,length(x),1);

D=zeros(N);
C = ones(1,N);  C(1) = 2;  C(N) = 2;
for i=1:N
    for j=1:N
        if (i~=j)
            D(i,j) = C(i)*(-1)^(i+j) / (C(j)*(x(i)-x(j)));
        elseif (i==1)
            D(i,i) = (2*(N-1)^2 + 1)/6;
        elseif (i==N)
            D(i,i) = -(2*(N-1)^2 + 1)/6;
        else
            D(i,i) = -x(i) / (2*(1-x(i)^2));
        end
    end
end