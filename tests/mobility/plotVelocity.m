clear all
close all
clc
fontSize=16;

material='W';
Tm=3695; % melting temperature[K]
mu=161.0e9; % shear modulus [Pa]
rho=19250.0;
cs=sqrt(mu/rho);
%materialFile='../../tutorials/DislocationDynamics/MaterialsLibrary/W.txt'
%materialFile='../../tutorials/DislocationDynamics/MaterialsLibrary/Cu.txt'
materialFile=['../../tutorials/DislocationDynamics/MaterialsLibrary/' material '.txt'];

system(['./mobility ' materialFile])

nT=101;

tSection=[300:100:800]/Tm;

dataS=load('velocityS.txt');
plotMobility(dataS,nT,fontSize,tSection,Tm,mu,cs)
print(gcf,[material '_prismScreMobility'], '-dpng', '-r300');

dataE=load('velocityE.txt');
plotMobility(dataE,nT,fontSize,tSection,Tm,mu,cs)
print(gcf,[material '_prismEdgewMobility'], '-dpng', '-r300');


function plotMobility(data,nT,fontSize,tSection,Tm,mu,cs)
datasize=size(data,1);

T=reshape(data(:,2),datasize/nT,nT);
S=reshape(data(:,1),datasize/nT,nT);
V=reshape(data(:,3),datasize/nT,nT);

figure
clf
hold on
surf(T,log10(S),log10(V+1e-10),'edgecolor','none','FaceAlpha',0.2)

isolevels=[1e-3];
for e=[-7:-1]
    for k=[10]
        isolevels=[isolevels k*10^e];
    end
end
%isolevels=[BrunnerV/cs isolevels];
%isolevels=[0.05:0.05:1]*cs
[C1,h1]=contour(T,log10(S),V,isolevels,'k','Linewidth',1);
clabel(C1,h1,'FontSize',fontSize,'Color','k','labelspacing', 700)
xlabel('T/Tm','FontSize',fontSize)
ylabel('log_{10}(\tau/\mu)','FontSize',fontSize)
grid on

% Const Temp sections are in columns of V
figure
hold on
tU=unique(data(:,2));
sU=unique(data(:,1));
for k=1:length(tSection)
tSec=tSection(k);
dt=abs(tU-tSec);
tID=find(dt==min(dt));
plot(sU*mu*1e-6,V(:,tID)*cs,'linewidth',2)
end
legend(num2str([tSection]'*Tm),'location','southeast')
xlabel('shear stress [MPa]')
ylabel('velocity [m/s]')
grid on
set(gca,'FontSize',16)
end


