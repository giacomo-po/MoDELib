close all
clear all
clc

fontsize=14;

Nx=201;
Ny=201;

Sn=load('S/S_0.txt');
Sa=load('S/S_1.txt');

X=reshape(Sn(:,1),Nx,Ny);
Y=reshape(Sn(:,2),Nx,Ny);
s11n=reshape(Sn(:,4),Nx,Ny);
s12n=reshape(Sn(:,5),Nx,Ny);
s22n=reshape(Sn(:,8),Nx,Ny);
s11a=reshape(Sa(:,4),Nx,Ny);
s12a=reshape(Sa(:,5),Nx,Ny);
s22a=reshape(Sa(:,8),Nx,Ny);

%% Surf Plots
figure(1)
clf
%plot3(X,Y,s11n)
surf(X,Y,s11n,'EdgeAlpha',0.05)

figure(2)
clf
surf(X,Y,s12n,'EdgeAlpha',0.05)

figure(3)
clf
surf(X,Y,s22n,'EdgeAlpha',0.05)

%% Plot on axes
figure(4)
clf
hold on
grid on
id=find(Y==0);
plot(X(id),s11n(id),'r')
plot(X(id),s11a(id),'r--')
plot(X(id),s22n(id),'b')
plot(X(id),s22a(id),'b--')
plot(X(id),s12n(id),'k')
plot(X(id),s12a(id),'k--')
axis([min(X(id)) max(X(id)) -0.2 0.2])
xlabel('x_1 / b','Fontsize',fontsize)
ylabel('stress / \mu','Fontsize',fontsize)
legend('\sigma_{11} (numerical)','\sigma_{11} (analytical)','\sigma_{22} (numerical)','\sigma_{22} (analytical)','\sigma_{12} (numerical)','\sigma_{12} (analytical)')
%print(gcf, '-depsc', 'sigma_edge_x1_Cai')

figure(5)
clf
hold on
grid on
id=find(X==0);
plot(Y(id),s11n(id),'r')
plot(Y(id),s11a(id),'r--')
plot(Y(id),s22n(id),'b')
plot(Y(id),s22a(id),'b--')
plot(Y(id),s12n(id),'k')
plot(Y(id),s12a(id),'k--')
axis([min(Y(id)) max(Y(id)) -0.2 0.2])
xlabel('x_2 / b','Fontsize',fontsize)
ylabel('stress / \mu','Fontsize',fontsize)
legend('\sigma_{11} (numerical)','\sigma_{11} (analytical)','\sigma_{22} (numerical)','\sigma_{22} (analytical)','\sigma_{12} (numerical)','\sigma_{12} (analytical)')
%print(gcf, '-depsc', 'sigma_edge_x2_Cai')
