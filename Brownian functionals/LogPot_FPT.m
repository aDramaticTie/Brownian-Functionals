%Distribution of the first passage time
%Brownian motion in log potential

tic
clear
clc
close all
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesFontSize',18);
pos=[0.11 0.11 0.85 0.85];

N=1e6; %number of paths
dt=1e-3;
NT=1e4; %number of time steps
T=dt*NT;
U=1.4;
D=1;
e=U/(2*D);
nu=0.5+e;
Init=1; %you can change the starting point
fpt=zeros(1,N);
for j=1:N
    %evolves the trajectory
    x=Init;
    tmp=0;
    fpt_r=0;
    for i=1:NT
        %We use the weak order 2 Runge-Kutta method
        dy=-U/x*dt+sqrt(2*D*dt)*randn;
        y=x+dy;
        dx=0.5*(-U/y-U/x)*dt+sqrt(2*D*dt)*randn;
        x=x+dx;
        tmp=tmp+dt;
        %---if you use
        %---x=x+sqrt(D*dt)*randn;
        %---then the theoretical pdf is
        %---1/sqrt(2*pi*D*T)*exp(-x.^2/(2*D*T))
        if x<0
            fpt_r=tmp;
            break
        end
    end
    fpt(j)=fpt_r;
end

fpt=fpt(fpt~=0);
%Build the histogram, using histcounts
createFit_scaled_Scott(fpt,1,1,'o','r');
hold on;

tt=logspace(-2,log10(T),10000);
td=Init^2/(4*D);
ftp_th=1/gamma(nu)*(td./tt.^(1+1/nu)).^nu.*exp(-td./tt);
plot(tt,ftp_th,'-k','linewidth',2);
toc