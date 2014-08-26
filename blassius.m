%blassius equation solver
function [u,difu]=blassius(y,ymax);
y=fliplr(y);
eta=0:0.001:75;					%Full solver 1E4 points
f0 = [0,0,0.3320573362151946];
[eta,f] = ode45(@BlasiusFunc,eta,f0);
u=f(:,2);
difu=f(:,3);
u=interp1(eta,u,y);				% u and u` gets interpolated to the y grid
difu=interp1(eta,difu,y);

u=fliplr(u);
difu=fliplr(difu);

