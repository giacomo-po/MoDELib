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
figure(3)
hold on

temp=zeros(500,16);
for k=[0:499]
fileID=k;
V=load(['V/V_' num2str(fileID) '.txt']);

temp(k+1,:)=V(:,8)';

end
plot(temp)