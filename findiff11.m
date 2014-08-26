function [x,D] = findiff11(xlower,xupper,N)

d = (xupper - xlower)/N;
x = xlower:d:xupper;
x = reshape(x,length(x),1);
% D2=zeros(N+1,N+1);
D = zeros(N+1,N+1);

for i=1:N+1
    for j=1:N+1
        if (i - j == 1)
            D(i,j) = 1/(2*d);
        elseif (i - j == -1)
            D(i,j) = -1/(2*d);
        elseif (i==j)
            D(i,j) = 0;
        else
            D(i,j) = 0;
        end
    end
end