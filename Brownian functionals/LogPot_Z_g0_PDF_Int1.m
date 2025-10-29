%PDF, Brownian motion in log potential in an interval (0,b)
%You need very small dt (1e-5), otherwise the numerical result is
%overestimated.
%Here I consider only g=0:if bb<0 we have a full distrib; if bb>0 we
%have a CONDITIONAL distrib; bb=0 is a critical case
%For the data set bb=0, we need a correcting factor for the normalization
%The CF1 is 0.9364, it comes from the integral of the PDF from 0 to 10
%(which is the max value considered for Z)
%The CF2 is 0.572786, it comes from the integral of the PDF from 0 to 0.2

tic
clear
clc
close all
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesFontSize',18);
pos=[0.12 0.12 0.85 0.83];

N=1e3; %number of paths
dt=1e-5;
NT=1e7; %number of time steps
T=dt*NT;
U=-1;
D=1;
e=U/(2*D);
bb=1+2*e;
x0=0.9; %you can change the starting point
b=1;

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
        if x>b
            res=tmp;
            break
        end
        if x<a
            break
        end
        tmp=tmp+dt*(x^(-2)+(x-dx)^(-2))/2;
    end
    Z(j)=res;
end
Z=Z(Z~=0);
Pr=numel(Z)/N;
zet=Z(Z<2);
Prob=numel(zet)/N;
[y,x]=createFit_FD(zet,'*','g');
figure(1)
semilogx(x,y/Pr,'*g','Markersize',7);
hold on;
zz=logspace(-3,1,10000);
yy=log(b/x0)./sqrt(4*pi*D*zz.^3).*exp(-((D*abs(bb)*zz-log(b/x0)).^2)./(4*D*zz));
semilogx(zz,yy/Prob,'-k','linewidth',2);
xlabel('$z$','interpreter','latex','fontsize',21.5);
ylabel('$p(z,x_0)$','interpreter','latex','fontsize',21.5);
set(gca,'Position',pos);

toc