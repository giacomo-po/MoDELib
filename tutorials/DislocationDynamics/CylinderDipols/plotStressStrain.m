clc
close all
clear all

fontSize=16;
S=load('S/S_0.txt');
%S=S(11:end,:);

B=1e-4;
mu=76e9;
%Sy=70e6/48e9;
Burgers=0.2489e-9; % [m]

R=2127/2;
H=6*R;
A=pi*R^2;
V=A*H;
eDot33=0.2*1.0e-9;



step=S(:,1);
%time=S(:,2);
dt=S(:,3);
time=cumsum(dt);
nLength=S(:,4);
s33=S(:,5);
eP11=S(:,6);
eP12=S(:,7);
eP13=S(:,8);
eP21=S(:,9);
eP22=S(:,10);
eP23=S(:,11);
eP31=S(:,12);
eP32=S(:,13);
eP33=S(:,14);

eE33=cumsum(eDot33-eP33).*dt;
figure(1)
plot(time*eDot33*100,s33*mu*1e-6);
grid on
set(gca,'FontSize',fontSize)
xlabel('strain [%]','FontSize',fontSize)
ylabel('stress [MPa]','FontSize',fontSize)

%return
%%
figure(2)
plot(eP33);
set(gca,'FontSize',fontSize)
ylabel('epsilon_P')
grid on

%%
figure(3)
plot(time*eDot33*100,nLength/V/Burgers^2);
set(gca,'FontSize',fontSize)
grid on


%%
return
uz=S(:,2)./S(:,1);
fz=S(:,3);

ez=-uz/H;
sz=-fz/A;

figure(1)
plot(ez*100,sz*mu*1e-6)
xlabel(texlabel('epsilon_z [%]'),'FontSize',fontSize)
ylabel(texlabel('sigma_z [MPa]'),'FontSize',fontSize)
grid on



return


dt=S(:,3);
de33=eDot33*dt;
e33=cumsum(de33);
s33=S(:,5);
nl=S(:,4);

ep11=cumsum(S(:,6).*dt);
ep12=cumsum(S(:,7).*dt);
ep13=cumsum(S(:,8).*dt);
ep21=cumsum(S(:,9).*dt);
ep22=cumsum(S(:,10).*dt);
ep23=cumsum(S(:,11).*dt);
ep31=cumsum(S(:,12).*dt);
ep32=cumsum(S(:,13).*dt);
ep33=cumsum(S(:,14).*dt);

figure(1)
plot(e33*100,s33,'Linewidth',2)
%plot(e33,'Linewidth',2)
hold on
plot(e33*100,e33*0+Sy,'r','Linewidth',1)
grid on
xlabel(texlabel('total strain epsilon_{33} [%]'),'fontSize',fontSize)
ylabel(texlabel('stress sigma_{33}/mu [-]'),'fontSize',fontSize)


figure(2)
%plot(e33,nl/V/Burgers^2,'Linewidth',2)
plot(nl/V/Burgers^2,'Linewidth',2)
grid on
xlabel(texlabel('stress_{33}/mu [-]'),'fontSize',fontSize)
ylabel(texlabel('dislocation density [1/m^2]'),'fontSize',fontSize)

figure(3)
%plot(e33*100,[ep11 ep12 ep13 ep21 ep22 ep23 ep31 ep32 ep33]*100,'Linewidth',1)
plot(e33*100,[ep11 ep12 ep13 ep21 ep22 ep23*0 ep31 ep32 ep33]*100,'Linewidth',1)
grid on
legend('Ep11','Ep12','Ep13','Ep21','Ep22','Ep23','Ep31','Ep32','Ep33')
xlabel(texlabel('total strain epsilon_{33} [%]'),'fontSize',fontSize)
ylabel(texlabel('plastic strain [%]'),'fontSize',fontSize)

figure(4)
%plot(e33*100,S(:,6:14),'Linewidth',1)
plot(e33*100,S(:,6:14),'Linewidth',1)
grid on
%xlabel(texlabel('total strain epsilon_{33} [%]'),'fontSize',fontSize)
ylabel(texlabel('plastic strain rate [?]'),'fontSize',fontSize)
legend('dotEp11','dotEp12','dotEp13','dotEp21','dotEp22','dotEp23','dotEp31','dotEp32','dotEp33')
set(gca,'fontSize',fontSize)


return
%%

for p=0:size(S,1)
    figure(5)
hold on
P=load(['P/P_' num2str(p) '.txt']);
plot(ones(size(P,1),1)*p,P(:,2),'x')
plot(p,mean(P(:,2)),'r.')
Vmax(p+1)=max(P(:,2));
if p==87
nMax=find(P(:,2)==max(P(:,2)))
nID=P(nMax,1)
V=load(['V/V_' num2str(p) '.txt']);
nRow=find(V(:,1)==nID)
figure(6)
axis equal
hold on
plot3(V(:,2),V(:,3),V(:,4),'bx')
plot3(V(nRow,2),V(nRow,3),V(nRow,4),'mo')
end
end

Vmax=Vmax(2:end);
plot(Vmax,'m')

return

%%
cs=1.888642561174662e-02;
figure(7)
subplot(2,1,1)
plot(Vmax)
hold on
plot(Vmax*0+cs*0.1,'r')
grid on
vcs=find(Vmax>cs*0.1);

%figure(8)
subplot(2,1,2)
dx=dt.*Vmax';
plot(dx)
hold on
plot(vcs,dx(vcs),'or')
grid on