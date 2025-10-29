%Average Z, Brownian motion in log potential in an interval (a,b)
%You need very small dt (1e-5), otherwise the numerical result is
%overestimated. You also need 1e4 paths. For a single starting point,
%the time is about 20s

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
U=2;
D=1;
e=U/(2*D);
bb=1+2*e;
g=3;
nu=bb/g;
a=0.5;
b=2.5;

%define a vector of starting points
ns=10;
x0=linspace(a,b,ns+2);
x0=x0(2:ns+1);
m=zeros(1,numel(x0));
st=m;

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
            if x<a || x>b
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

xx0=linspace(a,b,1000);
ZD=xx0.^g/(D*g^2);
%here I write g*nu instead of 1+2e
fa=b^(g*nu)-xx0.^(g*nu);
fb=xx0.^(g*nu)-a^(g*nu);
den=(b^(g*nu)-a^(g*nu))*xx0.^g;
mth=ZD/(nu-1).*(1-(fa*a^g+fb*b^g)./den);
if e==-0.5
    mth=ZD.*((b^g*log(xx0/a)+a^g*log(b./xx0))./(xx0.^g*log(b/a))-1);
end
if nu==1
    mth=ZD.*g.*(b^g*(xx0.^g-a^g).*log(b./xx0)-...
        a^g*(b^g-xx0.^g).*log(xx0/a))./(xx0.^g*(b^g-a^g));
end

figure(1)
plot(xx0,mth,'-b','linewidth',2);
hold on
plot(x0,m,'ok','Markersize',7,'Markerfacecolor','b');
xlabel('$x_0$','fontsize',21.5,'interpreter','latex');
ylabel('$\langle\mathcal{Z}\rangle$','fontsize',21.5,'interpreter','latex');
set(gca,'position',pos);
toc