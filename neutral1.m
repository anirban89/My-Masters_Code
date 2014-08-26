
 beta = -0.4;
 k = 0:0.1:3;

a=eye([1 161]);
b=fliplr(a);
c=zeros([1 161]);

%--------------------------JACOBIAN DEFINITION


% [yb,d]=chebdif(161,1);
[yb,d]=cheb1(161);
[dnew,y]=jd(yb,d);
% [U,difu]=blassius(y);
% diffu=dnew*difu;			
D=dnew;
% U=diag(U);diffu=diag(diffu);
% temp1=ones([61,241]);
 U1 = (sech(y-1)).^2 + (sech(y+1)).^2;
 U = U1./(max(U1));
 U = diag(U); 
% diagU = diag (U);
%-------------------------- NEUTRAL LOOP ITERATION

for n=1:31
% for m=1:36


% A=1/Re(m)*(d^4-2*k^2*d^2+k^4*eye(161))-1i*alp(n)*U*(d^2-k^2*eye(161)) + 1i*alp(n)*diffu;
% B=d^2-k^2*eye(161);

A = U*k(n)*(D^2 - k(n)^2*eye(161)) - (D^2 * U - beta)*k(n);
B = D^2 - k(n)^2*eye(161) ;

%-------------------------- BOUNDARY CONDITIONS
% A(1,:)=a;

A(1,:)=a;
A(161,:)=b;

% A(199,:)=d(161,:);

B(1,:)=c;
B(161,:)=c;
% B(2,:)=c;
% B(199,:)=c;

%A=A(2:199,2:199);
%B=B(2:199,2:199);

%%%%%%

% [vcap,lambda]=eig(A,B);lambda=-1/1i*diag(lambda);
% r=real(lambda);im=imag(lambda);
% r=r/k;im=im/k;

% temp1(m,n)=max(im);

[svect,sig] = eig(A,B);
sigr = real(sig); sigi = imag(sig);
% g=find(abs(sigi)>=100);
% sigi(g)=[];sigr(g)=[];
% g=find(abs(sigr)>=100);
% sigi(g)=[];sigr(g)=[];
% g=find(sigr>0.99 & sigi>0);
% sigi(g)=[];sigr(g)=[];
% gr(n) = max(max(sigi));
temp = -50000;
for in = 1:161
    if sigi(in,in) < 100 && sigi(in,in) > temp
        temp = sigi(in,in);
        tempr = sigr(in,in);
        maxi = in;
    end
end
gi(n) = temp;
gr(n) = tempr;

% n/241*100
end

% gr1=transpose(gr);
% % v=[0 0];
% contour(beta,k,gr1);
plot(k,gr)
figure, plot(k,gi)

% title('Neutral stability curve');
% xlabel('b');
% ylabel('k');
