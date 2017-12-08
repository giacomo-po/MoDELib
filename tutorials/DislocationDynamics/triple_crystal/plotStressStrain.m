close all
clear all
clc

fontSize=16;

F=load('F/F_0.txt');

e33=F(:,end-9);
s33=F(:,end);

figure(1)
plot([0;e33*100],[0;s33],'-')
grid on
xlabel('\epsilon_{33} [%]')
ylabel('\sigma_{33} / \mu ')
set(gca,'FontSize',fontSize)
print(gcf,'-depsc','stressStrainTripleCrystal')