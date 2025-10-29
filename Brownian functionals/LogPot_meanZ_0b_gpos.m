%Average Z, Brownian motion in log potential in an interval (0,b)
%You need very small dt (1e-5), otherwise the numerical result is
%overestimated. You also need 1e4 paths. For a single starting point,
%the time is about 20s
%Here I consider only g>0

tic
clear
clc
close all
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesFontSize',17);
pos=[0.12 0.12 0.85 0.83];

N=1e4; %number of paths
dt=1e-5;
NT=1e6; %number of time steps
T=dt*NT;
U=2;
D=1;
e=U/(2*D);
bb=1+2*e;
g=3;
b=2;

%define a vector of starting points
ns=10;
x0=linspace(0,b,ns+2);
x0=x0(2:ns+1);
m=zeros(1,numel(x0));
st=m;

Z=zeros(1,N);
for k=1:numel(x0)
    Z=zeros(1,N);
    for j=1:N
        %evolves the trajectory
        x=x0(k);
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
            if x<0 || x>b
                res=tmp;
                break
            end
            tmp=tmp+dt*(x-dx/2)^(g-2);
        end
        Z(j)=res;
    end
    Z=Z(Z~=0);
    m(k)=mean(Z);
    st(k)=std(Z);
end

xx0=linspace(0,b,1000);
ZD=xx0.^g/(D*g^2);
mth=ZD*g/(g-bb).*((b./xx0).^(g-max(0,bb))-1);
if bb==g
    mth=ZD*g.*log(b./xx0);
end

figure(1)
plot(xx0,mth,'-g','linewidth',2);
hold on
plot(x0,m,'sk','Markersize',7,'Markerfacecolor','g');
xlabel('$x_0$','fontsize',21.5,'interpreter','latex');
ylabel('$\langle\mathcal{Z}\rangle$','fontsize',21.5,'interpreter','latex');
set(gca,'position',pos);

toc