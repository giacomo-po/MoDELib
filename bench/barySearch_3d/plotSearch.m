close all
clear all
clc

dim=3;

meshID=1;
%x=rand(1,2)*30
%x=[-0.05 -0.15 0.0]
x=[-7.071067811865473e+02  7.071067811865477e+02  2.088793431625061e+03]
guess=[1318 2157 2160 3528];
%guess=[1318 2157 3528 5033];
system(['./search ' num2str(meshID) ' ' num2str(x(1)) ' ' num2str(x(2)) ' ' num2str(x(3)) ' ' num2str(guess(1)) ' ' num2str(guess(2)) ' ' num2str(guess(3)) ' ' num2str(guess(4))])

N=load(['N/N_' num2str(meshID) '.txt']);
P=load('P/P_0.txt');
S=load('S/S_0.txt');

figure(1)
hold on
grid on
axis equal


for n=2:size(S,1)
IDs=S(n,[4:7])+1;
V=N(IDs,[2:4])';
for i=1:4
for j=1:4
    plot3([V(1,i) V(1,j)],[V(2,i) V(2,j)],[V(3,i) V(3,j)],'k')
end
end
end


plot3(S(1,1),S(1,2),S(1,3),'ro')
plot3(S(2:end,1),S(2:end,2),S(2:end,3),'Linewidth',2)