%--------------------------- FLAT PLATE FLOW
clear all;
%--------------------------- BASIC DEFINITIONS
alp=0.01:0.001:0.25;
Re = 250:5:550;
val=0*Re;
beta = 0*alp;
a=eye([1 200]);
b=fliplr(a);
c=zeros([1 200]);

%--------------------------JACOBIAN DEFINITION


[yb,d]=chebdif(200,1);
[dnew,y]=jd(yb,d);
[U,difu]=blassius(y);
diffu=dnew*difu;			
d=dnew;
U=diag(U);diffu=diag(diffu);
temp1=ones([61,241]);

%-------------------------- NEUTRAL LOOP ITERATION

%n for alp, m for Re
for n=1:241
for m=1:61

k=(alp(n)^2+beta(n)^2)^0.5;
A=1/Re(m)*(d^4-2*k^2*d^2+k^4*eye(200))-1i*alp(n)*U*(d^2-k^2*eye(200)) + 1i*alp(n)*diffu;
B=d^2-k^2*eye(200);



%-------------------------- BOUNDARY CONDITIONS
A(1,:)=a;
A(200,:)=b;
A(2,:)=d(1,:);
A(199,:)=d(200,:);

B(1,:)=c;
B(200,:)=c;
B(2,:)=c;
B(199,:)=c;

%A=A(2:199,2:199);
%B=B(2:199,2:199);

%%%%%%

[vcap,lambda]=eig(A,B);lambda=-1/1i*diag(lambda);
r=real(lambda);im=imag(lambda);
r=r/k;im=im/k;
g=find(abs(im)>=100);
im(g)=[];r(g)=[];
g=find(abs(r)>=100);
im(g)=[];r(g)=[];
g=find(r>0.99 & im>0);
im(g)=[];r(g)=[];
temp1(m,n)=max(im);
end
n/241*100
end

temp1=transpose(temp1);
v=[0 0];
contour(Re,alp,temp1,v);

title('Neutral stability curve');
xlabel('Re');
ylabel('Alpha');
