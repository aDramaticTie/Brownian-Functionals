%Equivalence between Langevin and heterogeneous diffusion


tic
clear
clc
close all
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesFontSize',18);
pos=[0.12 0.12 0.85 0.83];

N=1e6; %number of paths
dt=1e-4;
NT=1e8; %number of time steps
T=dt*NT;
D=1;
bb=1.5;
g1=1;
g2=2;
g3=3;
x0=1; %you can change the starting point
Z1=zeros(1,N);
Z2=zeros(1,N);
Z3=zeros(1,N);

exp1=1-1/bb;
exp2=1-2/bb;
parfor j=1:N
    %evolves the trajectory
    x=x0;
    tmp1=0;
    tmp2=0;
    tmp3=0;
    res1=0;
    res2=0;
    res3=0;
    for i=1:NT
        %We use the weak order 2 Runge-Kutta method
        % Del=sqrt(2*D*dt)*x^(1-1/bb);
        % y_p=x+Del;
        % y_m=max(0,x-Del);
        % y=max(0,x+Del*randn);
        % dx=1/4*sqrt(2*D*dt)*(y_p^(1-1/bb)+y_m^(1-1/bb)+2*x^(1-1/bb))*randn...
        %     +1/4*sqrt(2*D*dt)*(y_p^(1-1/bb)-y_m^(1-1/bb))*(randn^2-1);
          % Y=.5*(randn+randn/sqrt(3));
          % dx=sqrt(2*D*dt)*x^(1-1/bb)*randn...
          %     +0.5*2*D*dt*(1-1/bb)*x^(1-2/bb)*(randn^2-1)...
          %     -0.5*(2*D*dt)^(1.5)*(1-1/bb)*(1/bb)*x^(1-3/bb)*(randn-Y);
        dx=sqrt(2*D*dt)*x^(exp1)*randn+D*(1-1/bb)*x^(exp2)*dt;
        x=x+dx;
        %---if you use
        %---x=x+sqrt(D*dt)*randn;
        %---then the theoretical pdf is
        %---1/sqrt(2*pi*D*T)*exp(-x.^2/(2*D*T))
        if x<0
            res1=tmp1;
            res2=tmp2;
            res3=tmp3;
            break
        end
        tmp1=tmp1+dt*(x^(g1-2)+(x-dx)^(g1-2))/2;
        tmp2=tmp2+dt*(x^(g2-2)+(x-dx)^(g2-2))/2;
        tmp3=tmp3+dt*(x^(g3-2)+(x-dx)^(g3-2))/2;
    end
    Z1(j)=res1;
    Z2(j)=res2;
    Z3(j)=res3;
end
g=g1;
Z=Z1(Z1~=0);
g_=(g-2)*bb+2;
nu=1/g_;
ZD=x0^(g_/bb)/((g_/bb)^2*D);
zet_=Z;
%zet_=Z/ZD;
zet=zet_(zet_<70);
%zet=zet_(zet_<90);
Prob=numel(zet)/N;
%createFit_scaled_Scott(zet,1,1,'o','b');
createFit_FD(zet,'o','b');
hold on;
zz=logspace(-2,2,10000);
yy=(ZD)^nu/gamma(nu)*(zz).^(-1-nu).*exp(-ZD./zz);
%yy=1/gamma(nu)*(zz).^(-1-nu).*exp(-1./zz);
plot(zz,yy/Prob,'-k','Linewidth',2);
xlabel('$z$','interpreter','latex','fontsize',21.5);
ylabel('$g(z,y_0)$','interpreter','latex','fontsize',21.5);
set(gca,'Position',pos);
toc