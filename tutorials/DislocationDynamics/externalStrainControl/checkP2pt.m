close all
clear all
clc

fontSize=16;
frameID=0;

P=load(['P/P_' num2str(frameID) '.txt']);
sizeP=size(P)
P2pt=load(['P2pt/P_' num2str(frameID) '.txt']);
sizeP2pt=size(P2pt)

normP=sqrt(sum((P(:,7:9).*P(:,7:9))'))';
dP=P-P2pt;
normdP=sqrt(sum((dP(:,7:9).*dP(:,7:9))'))';
id=find(normP>1e-5);
x=sort(normdP(id)./normP(id)*100);
ids=find(diff(x));
figure
plot(x(ids),ids/length(x),'o')
grid on
xlabel('stress error [%]','FontSize',fontSize)
ylabel('cumulative distribution','FontSize',fontSize)
set(gca,'FontSize',fontSize)