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

dt=1e-5;
NT=1e6; %number of time steps
T=dt*NT;
D=1;
bb=3/2;
U=D*(bb-1);
a=0.5;
b=2.5;

x0=1;
x=zeros(1,NT);
x(1)=x0;
for j=1:NT
    %evolves the trajectory
    %We use the weak order 2 Runge-Kutta method
    dy=-U/x(j)*dt+sqrt(2*D*dt)*randn;
    y=x(j)+dy;
    dx=0.5*(-U/y-U/x(j))*dt+sqrt(2*D*dt)*randn;
    x(j+1)=x(j)+dx;
    if x(j+1)<a || x(j+1)>b
        break
    end
end
x=x(x~=0);
F=x.^(-1);
tt=(0:numel(x)-1)*dt;
plot(tt,x,'-k','linewidth',.5);
hold on
plot(tt,F,'-k','linewidth',.5);
xlabel('$t$','interpreter','latex','fontsize',21.5);
ylabel('$F[x(t)]$','interpreter','latex','fontsize',21.5);
set(gca,'Position',pos);
toc