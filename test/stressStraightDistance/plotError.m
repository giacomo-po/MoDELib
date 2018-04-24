clear all
%close all
clc

data=load('output.txt');
c=data(:,1);
d=data(:,2);
L0=data(:,3);
L1=data(:,4);
uMin=data(:,5);
u1=data(:,6);
u2=data(:,7);
u3=data(:,8);

eMax1=data(:,9);
eMax2=data(:,10);
eMax3=data(:,11);
eMax3C=data(:,12);

%e1=data(:,5);
%eMax=data(:,6);
%eMaxC=data(:,7);
%eMax3=data(:,8);
%eMax3C=data(:,9);

%% plot
%figure(1)
%clf
%loglog(d,e,'.')
%loglog(d,e1,'rx')

%grid on

%x=d./(L1.*((L0./L1).^2+1));
%x=d./(0.01*L0+L1);
x=d./(0.00*L0+L1);

figure(1)
clf
hold on
scatter(x,eMax1,2,'filled','MarkerFaceAlpha',7/8,'MarkerFaceColor','b')
scatter(x,eMax2,2,'filled','MarkerFaceAlpha',7/8,'MarkerFaceColor','k')
scatter(x,eMax3,2,'filled','MarkerFaceAlpha',7/8,'MarkerFaceColor','r')
scatter(x,eMax3C,2,'filled','MarkerFaceAlpha',7/8,'MarkerFaceColor','g')
set(gca,'XScale','log')
set(gca,'YScale','log')
grid on
grid on
legend('1-pts interp','2-pts interp','3-pts interp','3-pts interp C')
set(gca,'FontSize',16)
axis(10.^[0 4 -4 0])
%figure(2)
%clf
%loglog(l./d,e1,'r.')
%grid on
x0=d./L0;
x1=d./L1;

figure(2)
clf
plot3(x0,x1,eMax1,'b.')
%hold on
%plot3(d./L0,d./L1,eMax2,'m.')

set(gca,'XScale','log')
set(gca,'YScale','log')
set(gca,'ZScale','log')
grid on
xlabel('d/L0','FontSize',16)
ylabel('d/L1','FontSize',16)

%%
figure(3)
yyaxis left
plot(eMax1,u1,'ko')
set(gca,'XScale','log')
ylabel('u1','FontSize',16)
yyaxis right
plot(eMax1,uMin,'bx')
grid on
ylabel('uMin','FontSize',16)
xlabel('eMax1','FontSize',16)

figure(4)
yyaxis left
plot(eMax2,u2,'ko')
set(gca,'XScale','log')
ylabel('u2','FontSize',16)
yyaxis right
plot(eMax2,uMin,'bx')
grid on
ylabel('uMin','FontSize',16)
xlabel('eMax2','FontSize',16)

figure(5)
yyaxis left
plot(eMax3,u3,'ko')
set(gca,'XScale','log')
ylabel('u3','FontSize',16)
yyaxis right
plot(eMax3,uMin,'bx')
grid on
ylabel('uMin','FontSize',16)
xlabel('eMax3','FontSize',16)
