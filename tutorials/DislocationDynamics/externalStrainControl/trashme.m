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
Ps=load('Ps/P_100.txt');
Pq=load('P/P_100.txt');
%Pc=load('Pc/P_0.txt');
%Pcs=load('Pcs/P_0.txt');
%id=128:128+64
id=1:size(Pq,1);
%e=sqrt(sum(((Ps(id,7:9)-P(id,7:9)).^2)'));
%q=sqrt(sum(((P(id,7:9)).^2)'));
%r=sqrt(sum(((Ps(id,7:9)).^2)'));
% qx=Pq(id,7);
% rx=Ps(id,7);
% qy=Pq(id,8);
% ry=Ps(id,8);
% qz=Pq(id,9);
% rz=Ps(id,9);

clrs={'b','r','k'}
figure(1)
clf
hold on

for c=7:9
%plot(e./r,'.')
plot(Pq(id,c),clrs{c-6})
plot(Ps(id,c),['--' clrs{c-6}])
%plot(Pc(id,c),['-.' clrs{c-6}])
%plot(Pcs(id,c),['.' clrs{c-6}])

end

id=find(diff(Ps(:,3))<0);
id1=[1;id+1];
id2=[id;size(Ps,1)];
%id1=find(diff(Ps(:,1)));
text(id1,id1*0+1e-3,num2str(Ps(id1,1)))
%id2=find(diff(Ps(:,2)));
text(id2,id2*0-1e-3,num2str(Ps(id2,2)))


figure(2)
clf
hold on
%plot(e./r,'.')
plot(Pq(:,[7:9])-Ps(:,[7:9]),'.')
%plot(qz,'.')
%legend('straight','quadrature')
%plot(s,'.')