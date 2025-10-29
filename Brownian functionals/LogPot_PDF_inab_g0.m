%Brownian motion in log potential in an interval (a,b)
%PDF for nu=1/2

tic
clear
clc
close all
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesFontSize',18);
pos=[0.12 0.12 0.85 0.83];

N=1e4; %number of paths
dt=1e-5;
NT=1e6; %number of time steps
T=dt*NT;
D=1;
bb=0;
U=D*(bb-1);
a=0.5;
b=2.5;

x0=1;
Z=zeros(1,N);
parfor j=1:N
%evolves the trajectory
    x=x0;
    tmp=0;
    res=0;
    for i=1:NT
        %We use the weak order 2 Runge-Kutta method
        dy=-U/x*dt+sqrt(2*D*dt)*randn;
        y=x+dy;
        dx=0.5*(-U/y-U/x)*dt+sqrt(2*D*dt)*randn;
        x=x+dx;
        if x<a || x>b
            res=tmp;
            break
        end
        tmp=tmp+dt*(x^(-2)+(x-dx)^(-2))/2;
    end
    Z(j)=res;
end
Z=Z(Z~=0);
Pr=numel(Z)/N;
createFit_scaled_Scott(Z,1,1,'o','b');
hold on

la=log(x0/a);
lb=log(b/x0);
L=log(b/a)/sqrt(D);
zz=logspace(-2,1,1000);
yy=zeros(1,numel(zz));
k_max=100;

for k=1:k_max
    tmp=(-1)^(k+1)*k*exp(-zz*(k*pi/L)^2)*...
        ((x0/b)^(bb/2)*sin(pi*k*la/L/sqrt(D))+(x0/a)^(bb/2)*sin(pi*k*lb/L/sqrt(D)));
    yy=yy+tmp;
end

plot(zz,yy*(2*pi/L^2).*exp(-zz*D*bb^2/4),'-k','linewidth',2);
toc