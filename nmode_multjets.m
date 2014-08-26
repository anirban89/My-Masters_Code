clear all

clc



kmax = 3;
kmin = 0;
nk = 60;
dk = (kmax - kmin)/nk;
k = kmin:dk:kmax;
N1 = 200;

blower = -2.5;
bupper = 1;
nbeta = 70;
dbeta = (bupper - blower)/nbeta ;
beta = blower:dbeta:bupper;

% L_k = 18;
% lk = 11:1:L_k - 1;

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
tempvect =zeros(N1+1,1) ;
maxi = 0;
gr = zeros(nbeta+1,nk+1);
ps = zeros(nbeta+1,nk+1);
eigfn = zeros(N1+1,nbeta+1,nk+1);
gr1 = zeros(nk+1,nbeta+1);
% mgr = zeros(L_k,1);
% si = zeros(N1+1,L_k);
%%%----------------------------------------------------------

[y,D2] = findiff1(-20,20,N1);


%  for o = 11:L_k 
 
    Lk = 3 ;
    
%%% BACKGROUND VELOCITY PROFILE-------------------------------

U1 = (sech(y-Lk)).^2 + (sech(y+Lk)).^2;
U = U1./(max(U1));
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
%         tempvect = svect(:,in);
        maxi = in;
    end
end
gr(m,n) = temp;
ps(m,n) = tempr;
% eigfn(:,m,n) = tempvect;

end
end


[an,bn] = find(gr == max(max(gr)));
% mgr(o) = max(max(gr));
% si(:,o) = real(eigfn(:,an,bn));
figure, plot(y,real(eigfn(:,an,bn)))
% title(sprintf('eigenfunction for most unstable mode ; L_k = %d',Lk))

gr1=transpose(gr);
figure, contour(beta,k,gr1);
xlabel('\beta')
ylabel('k')
title('Lk = 1')
% title(sprintf('U = sech^2(y-L_k) + sech^2(y+L_k) ; L_k = %d',Lk))
colorbar
% end 
% figure,plot(lk,mgr,'+-')
% xlabel('L_k')
% ylabel('growth rate (\sigma)')
