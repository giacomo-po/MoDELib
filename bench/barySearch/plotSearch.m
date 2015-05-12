close all
clear all
clc

dim=2;

meshID=3;
%x=rand(1,2)*30
x=[-0.05 -0.15]

system(['./search ' num2str(meshID) ' ' num2str(x(1)) ' ' num2str(x(2))])

P=load('P/P_0.txt');
S=load('S/S_0.txt');

figure(1)
hold on
grid on
axis equal
for t=1:dim:size(P,1)
fill(P(t,1:dim+1),P(t+1,1:dim+1),'g')
end

plot(S(1,1),S(1,2),'ro')
plot(S(2:end,1),S(2:end,2),'Linewidth',2)