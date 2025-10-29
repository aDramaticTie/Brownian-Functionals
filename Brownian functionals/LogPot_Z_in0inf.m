%First-passage functional, Brownian motion in log potential

tic
clear
clc
close all
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesFontSize',17);
pos=[0.12 0.12 0.85 0.83];

N=1e3; %number of paths
dt=1e-5;
NT=1e7; %number of time steps
T=dt*NT;
U=.5;
D=1;
bb=1+U/D;
g=1;
nu=bb/g;
x0=1; %you can change the starting point
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
ZD=x0^g/(g^2*D);
% createFit_scaled_Scott(Z,1,1,'o','b');
% hold on;
% zz=logspace(-2,1,10000);
% yy=1/gamma(nu)*(zz/ZD).^(-1-nu).*exp(-ZD./zz)/ZD;
% plot(zz,yy,'-k','Linewidth',2);
% xlabel('$z$','interpreter','latex','fontsize',21.5);
% ylabel('$p(z,x_0)$','interpreter','latex','fontsize',21.5);
% set(gca,'position',pos);
xi=Z/ZD;
createFit_scaled_Scott(xi,1,1,'^','g');
%createFit_nBin(xi,95,'s','r');
hold on;
zz=logspace(-2,1,10000);
yy=1/gamma(nu)*(zz).^(-1-nu).*exp(-1./zz);
plot(zz,yy/Prob,'-k','Linewidth',2);
xlabel('$\xi$','interpreter','latex','fontsize',21.5);
ylabel('$P(\xi)$','interpreter','latex','fontsize',21.5);
set(gca,'position',pos);

toc