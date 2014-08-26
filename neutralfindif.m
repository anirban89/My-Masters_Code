clear all

clc



kmax = 3;
kmin = 0;
nk = 60;
dk = (kmax - kmin)/nk;
k = kmin:dk:kmax;
N1 = 600;

blower = -2.5;
bupper = 1;
nbeta = 70;
dbeta = (bupper - blower)/nbeta ;
beta = blower:dbeta:bupper;


%%% DECLARE VARIABLES-------------------------------------

a=eye([1 N1+1]);
b=fliplr(a);
c=zeros([1 N1+1]);

A = zeros(N1+1,N1+1);
B = zeros(N1+1,N1+1);
sig = zeros(N1+1,N1+1);
svect = zeros(N1+1,N1+1);
sigr = zeros(N1+1,N1+1);
sigi = zeros(N1+1,N1+1);
temp = 0;
tempr = 0;
tempvect =zeros(N1+1) ;
maxi = 0;
gr = zeros(nbeta+1,nk+1);
ps = zeros(nbeta+1,nk+1);
eigfn = zeros(N1+1,nbeta+1,nk+1);
gr1 = zeros(nk+1,nbeta+1);


%%%----------------------------------------------------------

[y,D2] = findiff1(-20,20,N1);


% for Lk = 1:18; 
    
%%% BACKGROUND VELOCITY PROFILE-------------------------------

% U = (sech(y-Lk)+sech(y+Lk)).*tanh(y);
U = (sech(y));
diagU = diag(U); 

%%%-----------------------------------------------------------


for m = 1:nbeta +1
    
for n = 1:nk+1
   

A = diagU*k(n)*(D2 - k(n)^2*eye(N1+1)) - diag((D2 * U - beta(m))*k(n));
B = D2 - k(n)^2*eye(N1+1) ;

A(1,:)=a;
A(end,:)=b;


B(1,:)=c;
B(end,:)=c;

[svect,sig] = eig(A,B);
sigr = real(sig); sigi = imag(sig);

temp = -50000;
for in = 1:N1+1
    if sigi(in,in) < 100 && sigi(in,in) > temp
        temp = sigi(in,in);
        tempr = sigr(in,in);
        tempvect = svect(:,in);
        maxi = in;
    end
end
gr(m,n) = temp;
ps(m,n) = tempr;
eigfn(:,m,n) = tempvect;

end
end


[an,bn] = find(gr == max(max(gr)));
figure, plot(y,real(eigfn(:,an,bn)))
title(sprintf('eigenfunction for most unstable mode ; L_k = %d',Lk))

gr1=transpose(gr);
figure, contour(beta,k,gr1);
xlabel('\beta')
ylabel('k')
title(sprintf('U = (sech(y-L_k) + sech(y+L_k))*tanh(y) ; L_k = %d',Lk))
colorbar
% end 
