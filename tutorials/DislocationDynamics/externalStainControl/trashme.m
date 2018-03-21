close all
clear all
clc

% V=load('V/V_0.txt');
% Vs=load('Vs/V_0.txt');
% 
% figure
% hold on
% vs=sqrt(sum((Vs(:,5:7).^2)'))
% v=sqrt(sum((V(:,5:7).^2)'))
% 
% plot(Vs(:,1),vs)
% plot(V(:,1),v)
% 
 legend('straight','quadrature')

%return
Ps=load('Ps/P_10.txt');
Pq=load('P/P_10.txt');
%id=128:128+64
id=1:size(Pq,1)
%e=sqrt(sum(((Ps(id,7:9)-P(id,7:9)).^2)'));
%q=sqrt(sum(((P(id,7:9)).^2)'));
%r=sqrt(sum(((Ps(id,7:9)).^2)'));
qx=Pq(id,7);
rx=Ps(id,7);
qy=Pq(id,8);
ry=Ps(id,8);
qz=Pq(id,9);
rz=Ps(id,9);

clrs={'b','r','k'}
figure(1)
clf
hold on

for c=7:9
%plot(e./r,'.')
plot(Pq(id,c),clrs{c-6})
plot(Ps(id,c),['--' clrs{c-6}])
end


figure(2)
clf
hold on
%plot(e./r,'.')
plot(Pq(id,[7:9])-Ps(id,[7:9]),'.')
%plot(qz,'.')
%legend('straight','quadrature')
%plot(s,'.')