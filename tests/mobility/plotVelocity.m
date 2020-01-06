clear all
close all
clc
fontSize=16;

materialFile='../../tutorials/DislocationDynamics/MaterialsLibrary/W.txt'
%materialFile='../../tutorials/DislocationDynamics/MaterialsLibrary/Cu.txt'

system(['./bccMobility ' materialFile])

data=load('velocity.txt');
nT=101;
datasize=size(data,1);

T=reshape(data(:,2),datasize/nT,nT);
S=reshape(data(:,1),datasize/nT,nT);
V=reshape(data(:,3),datasize/nT,nT);

figure(1)
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