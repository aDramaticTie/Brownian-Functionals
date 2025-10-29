%First-passage functional, Brownian motion in log potential

tic
clear
clc
close all
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesFontSize',18);
pos=[0.11 0.11 0.85 0.85];

N=1e2; %number of paths
dt=1e-5;
NT=1e7; %number of time steps
T=dt*NT;
U=1.8;
D=1;
bb=1+U/D;
g=3;
nu=(1+2*e)/g;
Init=1; %you can change the starting point
Z=zeros(1,N);
parfor j=1:N
    %evolves the trajectory
    x=Init;
    tmp=0;
    res=0;
    for i=1:NT
        %We use the weak order 2 Runge-Kutta method
        dy=-U/x*dt+sqrt(2*D*dt)*randn;
        y=x+dy;
        dx=0.5*(-U/y-U/x)*dt+sqrt(2*D*dt)*randn;
        x=x+dx;
        %---if you use
        %---x=x+sqrt(D*dt)*randn;
        %---then the theoretical pdf is
        %---1/sqrt(2*pi*D*T)*exp(-x.^2/(2*D*T))
        if x<0
            res=tmp;
            break
        end
        tmp=tmp+dt*(x^(g-2)+(x-dx)^(g-2))/2;
    end
    Z(j)=res;
end
Z=Z(Z~=0);
Prob=numel(Z)/N;

figure(1)
ZD=Init^g/(g^2*D);
xx=logspace(-3,2,10000);
yy=1/gamma(nu)*(ZD./xx.^(1+1/nu)).^nu.*exp(-ZD./xx);
plot(xx,yy,'-k','Linewidth',2);
hold on;
createFit_scaled_Scott(Z,1,1,'s','r');

toc