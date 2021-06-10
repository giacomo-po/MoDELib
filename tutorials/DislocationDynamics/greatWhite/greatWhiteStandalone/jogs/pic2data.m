close all
clear all
clc

[filename, pathname] = uigetfile('*.jpg', 'Pick a jpg image');

%% read and plot image
I=imread([pathname filename]);

figure(1)
image(I)

%% Get axis
A=ginput(3);
dx_p=A(2,1)-A(1,1);
dy_p=A(3,2)-A(1,2);
x0=0;
y0=0;
dx_r=70e-9;
dy_r=0.25;

%% dataset 1
figure(1)
d1=ginput();


d1_0=d1-repmat(A(1,:),size(d1,1),1);
d1_real=[d1_0(:,1)/dx_p*dx_r+x0 d1_0(:,2)/dy_p*dy_r+y0]

%% plot
figure(2)
plot(d1_real(:,1),d1_real(:,2),'o','Linewidth',2)
grid on
axis([x0 max(d1_real(:,1)) 0 max(d1_real(:,2)) ])


