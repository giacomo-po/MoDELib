clc
clear all
close all

S=load('D/D_0.txt');

%% Plot function
comp={'x_1','value'};
figure(2)
hold on

clrCol=2;
sMax=max(S(:,clrCol))
sMin=min(S(:,clrCol))
clrID=round((S(:,clrCol)-sMin)/(sMax-sMin)*size(colormap,1));

for f=1:2:size(S,1)
    plot(S(f:f+1,1),S(f:f+1,2))
end
grid on
xlabel('X_1','FontSize',14)
ylabel('value','FontSize',14)
