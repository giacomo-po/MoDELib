% N=load('N/N_3.txt');
% T=load('T/T_3.txt');
% 
% T1=[T zeros(size(T,1),1)];
% for t=1:size(T,1)
%     vIDs=T(t,2:end)+1;
% if (mean(N(vIDs,3))<-300)
%         T1(t,5)=0;
% elseif (mean(N(vIDs,3))>=-300 && mean(N(vIDs,2)>150))
%         T1(t,5)=1;
%     else
%         T1(t,5)=2;
%     end
% end
% dlmwrite('myFile.txt', T1,'delimiter','\t');
%return

close all
clear all
clc

dim=2;

meshID=3;
%x=[-0.05 -0.15]
%x=[200 -600]
x=[1.493572e+02 -3.036185e+02]
regionID=2;

system(['./search ' num2str(meshID) ' ' num2str(x(1)) ' ' num2str(x(2)) ' ' num2str(regionID)])

P=load('P/P_0.txt');
S=load('S/S_0.txt');

figure(1)
hold on
grid on
axis equal
clrs={'r','g','b'}
for t=1:dim:size(P,1)
fill(P(t,1:dim+1),P(t+1,1:dim+1),clrs{P(t,dim+2)+1})
end

plot(S(1,1),S(1,2),'mo','Linewidth',2)
plot(S(2:end,1),S(2:end,2),'w','Linewidth',2)