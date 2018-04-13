%close all
clear all
clc


data=[1 0 2.5e-05 0 1e-06
10 0 7.5e-05 0 4e-06
100 0 0.000122 0 0.000372
1000 5 0.001503 5 0.02618
10000 619 0.030126 621 2.48689
100000 61291 2.0629 61412 251.439];

figure
clf 
hold on
plot(data(:,1),data(:,3),'o-','Linewidth',2)
plot(data(:,1),data(:,5),'o-','Linewidth',2)
set(gca,'XScale','log')
set(gca,'YScale','log')
xlabel('number of lines')
ylabel('cpu time [sec]')
grid on
legend('sweep-line','pair-intersection')