clear all
clc


kmax = 3;
kmin = 0;
nk = 60;
dk = (kmax - kmin)/nk;
k = kmin:dk:kmax;
Ny = 500;
Nx = 100;

nd = 5;        % index of the maximum separation of jets

blower = -2;
bupper = 2;
nbeta = 80;
dbeta = (bupper - blower)/nbeta ;
beta = blower:dbeta:bupper;

%%-------------------------------------
%% INITIALIZATION OF VALUES
%--------------------------------------
a=eye([1 Ny+1]);
b=fliplr(a);
c=zeros([1 Ny+1]);

A = zeros(Ny+1,Ny+1);
B = zeros(Ny+1,Ny+1);
sig = zeros(Ny+1,Ny+1);
svect = zeros(Ny+1,Ny+1);
sigr = zeros(Ny+1,Ny+1);
sigi = zeros(Ny+1,Ny+1);
temp = 0;
tempr = 0;
tempvect =zeros(Ny+1,1) ;
maxi = 0;
gr = zeros(nbeta+1,nk+1);
ps = zeros(nbeta+1,nk+1);
eigfn = zeros(Ny+1,nbeta+1,nk+1);
gr1 = zeros(nk+1,nbeta+1);
mgr = zeros(nd,1);
% si = zeros(Ny+1,L_k);
si1 = zeros(Ny+1,1);

%%---------------------------------------


[y,Dy,D2y] = findiff12(-20,20,Ny);
% x = -pi:0.01:pi;
% x = transpose(x);
[x,Dx] = findiff11(-pi,pi,Nx);

z1 = sin(x);

for o = 1:nd;

 Lk = o - 1;               % Lk = actual distance between two jets

%%-----------------------------------------
%% BACKGROUND VELOCITY PROFILE
%%-----------------------------------------
 
U1 =  (sech(y-Lk) + sech(y+Lk)).*tanh(y);
U = U1/(max(U1));;
diagU = diag(U);


%%-------------------------------------------------------------------------


for m = 1:nbeta +1
    
for n = 1:nk+1
   

A = diagU*k(n)*(D2y - k(n)^2*eye(Ny+1)) - diag((D2y * U - beta(m))*k(n));
B = D2y - k(n)^2*eye(Ny+1) ;

A(1,:)=a;
A(end,:)=b;


B(1,:)=c;
B(end,:)=c;

[svect,sig] = eig(A,B);
sigr = real(sig); sigi = imag(sig);

temp = -50000;
for in = 1:Ny+1
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
mgr(o) = max(max(gr));
%si(:,o) = real(eigfn(:,an,bn));
si1 = real(eigfn(:,an,bn));
figure, plot(y,real(eigfn(:,an,bn)))
title(sprintf('eigenfunction for most unstable mode U = (sech(y-L_k) + sech(y+L_k))tanh(y) ; L_k = %d',Lk))
%title('eigenfunction for most unstable mode; L_k = 3')

gr1=transpose(gr);
figure, contour(beta,k,gr1);
xlabel('\beta')
ylabel('k')
title('Lk = 1')
title(sprintf('U = (sech(y-L_k) + sech(y+L_k))tanh(y) ; L_k = %d',Lk))
% title('U = sech^2(y-L_k) + sech^2(y+L_k) ; L_k = 3')
colorbar

psi = (real(eigfn(:,an,bn))) * transpose(z1);
u_1 = - Dy * psi;
v_1 =  transpose(Dx * transpose(psi));

figure,contourf(x,y,psi)
colorbar
title(sprintf('streamfunction for two jets with basic velocity profile: U = (sech(y-L_k) + sech(y+L_k))tanh(y) ; L_k = %d',Lk))

X = x(1:5:end);
Y = y(1:5:end);
u = u_1(1:5:end,1:5:end);
v = v_1(1:5:end,1:5:end);


 figure,quiver(X,Y,u,v)
%  axis([-3.2 3.2 -(Lk+7) (Lk+7)])
 title(sprintf('velocity field for basic velocity profile U = (sech(y-L_k) + sech(y+L_k))tanh(y); L_k = %d ' , Lk))

end 
 figure,plot(Lk,mgr,'+-')
 xlabel('L_k')
 ylabel('growth rate (\sigma)')
