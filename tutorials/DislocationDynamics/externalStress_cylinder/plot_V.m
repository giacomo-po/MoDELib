close all
clear all
clc

figure(1)
hold on

for k=[50]
fileID=k;
V=load(['V/V_' num2str(fileID) '.txt']);

v=V(:,[5:7]);
vN=sqrt(sum((v.*v)'))';
vrf=V(:,8);

[vN,I] = sort(vN);
vrf=vrf(I);

plot(vN,vrf,'Color',rand(1,3))
end

%% 

nNodes=8;
nFiles=200;
vFilter=zeros(nFiles,nNodes);
v=zeros(nFiles,nNodes);
for k=[0:nFiles-1]
fileID=k;
V=load(['V/V_' num2str(fileID) '.txt']);

vFilter(k+1,:)=V(:,8)';
v(k+1,:)=sqrt(sum((V(:,[5:7]).*V(:,[5:7]))'))';
end

activeID=[0 2 4 6] ;
%activeID=[0:15];
%activeID=7;

figure(3)
hold on
plot(vFilter(:,activeID+1))
legend(num2str(activeID'))
grid on

figure(4)
hold on
plot(v(:,activeID+1))
legend(num2str(activeID'))
grid on
