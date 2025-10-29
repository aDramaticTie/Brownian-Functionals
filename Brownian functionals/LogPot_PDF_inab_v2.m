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
g=-3;
bb=-g/2;
U=D*(bb-1);
nu=bb/g;
a=0.5;
b=2.5;

x0=1;
Z=zeros(1,N);
Z0=zeros(1,N);
parfor j=1:N
%evolves the trajectory
    x=x0;
    tmp=0;
    tmp0=0;
    res=0;
    res0=0;
    for i=1:NT
        %We use the weak order 2 Runge-Kutta method
        dy=-U/x*dt+sqrt(2*D*dt)*randn;
        y=x+dy;
        dx=0.5*(-U/y-U/x)*dt+sqrt(2*D*dt)*randn;
        x=x+dx;
        if x<a || x>b
            res=tmp;
            res0=tmp0;
            break
        end
        tmp=tmp+dt*(x^(g-2)+(x-dx)^(g-2))/2;
        tmp0=tmp0+dt*(x^(-2)+(x-dx)^(-2))/2;
    end
    Z(j)=res;
    Z0(j)=res0;
end
zet=Z0;
if zet==Z0
    c=0;
else
    c=1;
end
zet=zet(zet<2);
Pr=numel(zet)/N;
createFit_scaled_Scott(zet,1,1,'^','g');
hold on
zz=logspace(-3,1,1000);
yy=zeros(1,numel(zz));
k_max=100;
if c==1
    L=2*sqrt((a^(g/2)-b^(g/2))^2/(D*g^2));
    la=2*sqrt((x0^(g/2)-a^(g/2))^2/(D*g^2));
    lb=2*sqrt((x0^(g/2)-b^(g/2))^2/(D*g^2));
    if nu>0
        for k=1:k_max
            tmp=(-1)^(k+1)*k*exp(-zz*(k*pi/L)^2)*(sin(pi*k*la/L)+sin(pi*k*lb/L));
            yy=yy+tmp;
        end
    else
        for k=1:k_max
            tmp=(-1)^(k+1)*k*exp(-zz*(k*pi/L)^2)*((x0/b)^bb*sin(pi*k*la/L)...
                +(x0/a)^bb*sin(pi*k*lb/L));
            yy=yy+tmp;
        end
    end
    plot(zz,yy*(2*pi/L^2),'-k','linewidth',2);
else
    la=log(x0/a);
    lb=log(b/x0);
    L=log(b/a)/sqrt(D);
    for k=1:k_max
        tmp=(-1)^(k+1)*k*exp(-zz*(k*pi/L)^2)*...
        ((x0/b)^(bb/2)*sin(pi*k*la/L/sqrt(D))+(x0/a)^(bb/2)*sin(pi*k*lb/L/sqrt(D)));
        yy=yy+tmp;
    end
    plot(zz,yy*(2*pi/L^2).*exp(-zz*D*bb^2/4),'-k','linewidth',2);
end
xlabel('$z$','interpreter','latex','fontsize',21.5);
ylabel('$p(z,x_0)$','interpreter','latex','fontsize',21.5);
set(gca,'position',pos);
toc