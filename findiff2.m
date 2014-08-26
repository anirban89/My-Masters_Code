function [x,D2] = findiff2(xlower,xupper,N)

d = (xupper - xlower)/N;
x = xlower:d:xupper;
x = reshape(x,length(x),1);
D2=zeros(N);

for i=1:N+1
    for j=1:N+1
        if (abs(i - j) == 1)
            D2(i,j) = 1/(5*d^2);
        elseif (abs(i - j) == 2)
            D2(i,j) = 1/(5*d^2);
        elseif (i==j)
            D2(i,j) = -4/(5*d^2);
        else
            D2(i,j) = 0;
        end
    end
end