%Average Z, Brownian motion in log potential in an interval (0,b)
%You need very small dt (1e-5), otherwise the numerical result is
%overestimated. You also need 1e4 paths. For a single starting point,
%the time is about 20s
%Here I consider only g<0:if bb<g we have a finite mean value; if bb>-g we
%have a finite CONDITIONAL mean value; in between we have an infinite one.

%pos_i=[.5 .5 .4 .4]; %internal position for insets
%commands for insets: ax2=axes('Position',pos_i); creates inset
%now you can plot in the inset. Remember hold(ax2,'on');

tic
clear
clc
close all
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesFontSize',17);
pos=[0.12 0.12 0.85 0.83];
pos_i=[.5 .5 .4 .4];

N=1e5; %number of paths
dt=1e-5;
NT=1e7; %number of time steps
T=dt*NT;
U=-1.2;
D=1;
e=U/(2*D);
bb=1+2*e;
g=-.5;
b=2;

%define a vector of starting points
ns=10;
x0=linspace(0,b,ns+2);
x0=x0(2:ns+1);
m=zeros(1,numel(x0));
ste=m;  %standard error: when you estimate the mean, you consider the st
        %error, which is the st deviation divided by 
Eb=m;   %Pr of leaving from b

Z=zeros(1,N);

%First do the case with finite mean value
if bb<g || bb>-g    
    for k=1:numel(x0)
        Z=zeros(1,N);
        x=x0(k);
        parfor j=1:N
            %evolves the trajectory
            x_=x;
            tmp=0;
            res=0;
            for i=1:NT
                %We use the weak order 2 Runge-Kutta method
                dy=-U/x_*dt+sqrt(2*D*dt)*randn;
                y=x_+dy;
                dx=0.5*(-U/y-U/x_)*dt+sqrt(2*D*dt)*randn;
                x_=x_+dx;
                %---if you use
                %---x=x+sqrt(D*dt)*randn;
                %---then the theoretical pdf is
                %---1/sqrt(2*pi*D*T)*exp(-x.^2/(2*D*T))
                if x_>b
                    res=tmp;
                    break
                end
                if x_<0
                    break
                end
                %tmp=tmp+dt*(x_-dx/2)^(g-2);
                %metodo trapezi
                tmp=tmp+dt*(x_^(g-2)+(x_-dx)^(g-2))/2;
            end
            Z(j)=res;
        end
        Z=Z(Z~=0);
        Eb(k)=numel(Z)/N;
        m(k)=mean(Z);
        ste(k)=std(Z)/sqrt(numel(Z));
    end
    xx0=linspace(0,b,1000);
    ZD=xx0.^g/(D*g^2);
    mth=ZD*g/(g+abs(bb)).*((b./xx0).^g-1);
    Ebth=(xx0/b).^bb;

figure(1)
semilogy(xx0,mth,'-r','linewidth',2);
hold on
semilogy(x0,m,'sk','Markersize',7,'Markerfacecolor','r');
xlabel('$x_0$','fontsize',21.5,'interpreter','latex');
ylabel('$\langle\mathcal{Z}_b\rangle$','fontsize',21.5,'interpreter','latex');
set(gca,'position',pos);

%Only useful when bb>0
figure(2)
plot(xx0,Ebth,'-k','linewidth',2);
hold on
plot(x0,Eb,'dk','Markersize',7,'Markerfacecolor','r');
xlabel('$x_0$','fontsize',21.5,'interpreter','latex');
ylabel('$\mathcal{E}_b$','fontsize',21.5,'interpreter','latex');
set(gca,'position',pos);

end
if g<bb && bb<-g
    Z=zeros(1,N);
    x0=0.5;
    parfor j=1:N
        %evolves the trajectory
        x_=x0;
        tmp=0;
        res=0;
        for i=1:NT
            %We use the weak order 2 Runge-Kutta method
            dy=-U/x_*dt+sqrt(2*D*dt)*randn;
            y=x_+dy;
            dx=0.5*(-U/y-U/x_)*dt+sqrt(2*D*dt)*randn;
            x_=x_+dx;
            %---if you use
            %---x=x+sqrt(D*dt)*randn;
            %---then the theoretical pdf is
            %---1/sqrt(2*pi*D*T)*exp(-x.^2/(2*D*T))
            if x_>b
                res=tmp;
                break
            end
            if x_<0
                break
            end
            %tmp=tmp+dt*(x_-dx/2)^(g-2);
            %metodo trapezi
            tmp=tmp+dt*(x_^(g-2)+(x_-dx)^(g-2))/2;
        end
        Z(j)=res;
    end
    Z=Z(Z~=0);
end

toc