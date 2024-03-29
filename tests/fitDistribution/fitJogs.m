clear all
close all
clc
fontSize=16;


%% Jog Length
jogL_644K=[0.000000125   0.045774647887324
   0.000000150   0.136443661971831
   0.000000175   0.090669014084507
   0.000000200  0
   0.000000225000000   0.227112676056338
   0.000000250000000   0.091549295774648
   0.000000275   0.136443661971831
   0.000000300000000   0.090669014084507
   0.000000325000000   0.090669014084507
   0.000000350000000   0.045774647887324
      0.000000375  0
   0.000000400  0
   0.000000425000000   0.045774647887324];
jogL_644K(:,2)=jogL_644K(:,2)/sum(jogL_644K(:,2)); % renormalize fractions

xMax=500e-9;
dx=25e-9;
x=[0:dx:xMax];

% Normal distribution
m=238e-9; %initial guess for mean
s=50e-9; %initial guess for std
f=@(p,x) 1./(sqrt(2*pi)*p(2)).*exp(-(x-p(1)).^2/(2*p(2)^2)); % p(1)=mean, p(2)=std
p0=[m,s];
[p,resnorm,~,exitflag,output] = lsqcurvefit(f,p0,jogL_644K(:,1),jogL_644K(:,2)/dx)

% Lognormal distribution
m=log(238e-9); %initial guess for mean
s=1; %initial guess for std
g=@(p,x) 1./(sqrt(2*pi)*p(2)*x).*exp(-(log(x)-p(1)).^2/(2*p(2)^2)); % p(1)=mean, p(2)=std
q0=[m,s];
[q,resnorm,~,exitflag,output] = lsqcurvefit(g,q0,jogL_644K(:,1),jogL_644K(:,2)/dx)


figure(1)
bar(jogL_644K(:,1),jogL_644K(:,2),1)
hold on
plot(x,f(p,x)*dx,'Linewidth',2)
plot(x,g(q,x)*dx,'Linewidth',2)
grid on
ylabel('fraction')
xlabel('jog length [m]')
set(gca,'FontSize',fontSize)
legend('experimental',['normal: mean=' num2str(p(1)*10^9) 'nm, std=' num2str(p(2)*10^9) 'nm'],['log-normal: m=' num2str(q(1)) ', s=' num2str(q(2)) ])

%% Jog heigth
jogH_644K=[
   0.000000017   0.125461254612546
   0.000000021   0.125461254612546
   0.000000025   0.250922509225092
   0.000000029   0.126383763837638
   0.000000033   0.188191881918819
   0.000000037   0.124538745387454
   0.000000053   0.062730627306273];
jogH_644K(:,2)=jogH_644K(:,2)/sum(jogH_644K(:,2)); % renormalize fractions
xMax=70e-9;
dx=4e-9;
x=[1e-9:dx:xMax];

m=27e-9; %initial guess for mean
s=5e-9; %initial guess for std
f=@(p,x) 1./(sqrt(2*pi)*p(2)).*exp(-(x-p(1)).^2/(2*p(2)^2)); % p(1)=mean, p(2)=std
p0=[m,s];
[p,resnorm,~,exitflag,output] = lsqcurvefit(f,p0,jogH_644K(:,1),jogH_644K(:,2)/dx)

% Lognormal distribution
m=log(27e-9); %initial guess for mean
s=0.5; %initial guess for std
g=@(p,x) 1./(sqrt(2*pi)*p(2)*x).*exp(-(log(x)-p(1)).^2/(2*p(2)^2)); % p(1)=mean, p(2)=std
q0=[m,s];
[q,resnorm,~,exitflag,output] = lsqcurvefit(g,q0,jogH_644K(:,1),jogH_644K(:,2)/dx)

figure(2)
bar(jogH_644K(:,1),jogH_644K(:,2),1)
hold on
plot(x,f(p,x)*dx,'Linewidth',2)
plot(x,g(q,x)*dx,'Linewidth',2)
grid on
ylabel('fraction')
xlabel('jog height [m]')
set(gca,'FontSize',fontSize)
legend('experimental',['normal: mean=' num2str(p(1)*10^9) 'nm, std=' num2str(p(2)*10^9) 'nm'],['log-normal: m=' num2str(q(1)) ', s=' num2str(q(2)) ])