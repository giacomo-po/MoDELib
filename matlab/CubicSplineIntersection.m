function CubicSplineIntersection(B1,B2)




%B2=[472.8 493.1 519.5 543.6;
%1.426 1.215 3.315 -13.94];

%B1=rand(3,4);
%B2=rand(3,4);


%Zx=ones(3,4);
%Zx(1,:)=0;

%Zy=ones(3,4);
%Zy(2,:)=0;


figure(1)
plotSpline(B1,'b')
plotSpline(B2,'r')
%plotSpline(B1.*Zx,'b')
%plotSpline(B2.*Zx,'r')
%plotSpline(B1.*Zy,'b')
%plotSpline(B2.*Zy,'r')


[sx,tx]=PlanarIntersection(B1([1 2],:),B2([1 2],:))
R1=position(B1,tx);
R2=position(B2,sx);

figure(1)
plot(R1(1),R1(2),'go','Linewidth',2)
figure(1)
plot(R2(1),R2(2),'ms','Linewidth',2)
%[sy,ty]=PlanarIntersection(B1([1 3],:),B2([1 3],:))


%roots = intersect ([sx' tx'],[sy' ty'],'rows' )

function [s,t]=PlanarIntersection(B1,B2)



[s,t]=findIntersections(B1,B2);
%s=zeros(1,9);
%t=zeros(1,9);

epsilon=0;
commonroots=(abs(s-0.5)<0.5+epsilon).*(abs(t-0.5)<0.5+epsilon);

s=s(find(commonroots));
t=t(find(commonroots));

%s(2)=0.5757;
%t(2)=0.5757;

%s=s(find(real(s)>=0 ))
%s=s(find(real(s)<=1 ))
%s=s(find(imag(s)==0 ))
%s(1)
%s(2)
%s(4)
%s=s(find(real(s)<=1 ))
%real(s)
%R1=position(B1,real(t));
%R2=position(B2,real(s));

R1=real(position(B1,t));
R2=real(position(B2,s));

% figure(1)
% %clf
% %subplot(1,2,1)
% plotSpline(B1,'r')
% plotSpline(B2,'m')
% %v=axis;
% plot(R1(1,:),R1(2,:),'bx','Linewidth',2)
% plot(R2(1,:),R2(2,:),'m.','Linewidth',2)
% plot([R1(1,:); R2(1,:)],[R1(2,:); R2(2,:)],'g','Linewidth',1)
% %axis(v)
% %plot()

figure(2)
subplot(1,2,2)
np=100;
theta=[0:np]/np*2*pi;
cx=0.5+(0.5+epsilon)*cos(theta);
cy=(0.5+epsilon)*sin(theta);
d=0.2;
plot([0 1],[0 0],'g','Linewidth',2)
hold on
patch(cx,cy,'y','FaceAlpha',0.3)
plot(real(s),imag(s),'bx','Linewidth',2)
plot(real(t),imag(t),'m.','Linewidth',2)
axis([-d 1+d -0.5-d 0.5+d])
grid on
axis equal

%end

function [s,t]=findIntersections(B1,B2)
B1E=[B1;ones(1,4)];
B2E=[B2;ones(1,4)];

M0=SubMatrixPolynomialFunction(B1E, B2E, 0);
M1=SubMatrixPolynomialFunction(B1E, B2E, 1);
M2=SubMatrixPolynomialFunction(B1E, B2E, 2);
M3=SubMatrixPolynomialFunction(B1E, B2E, 3);
%return
%[V,D] = eigenSolver(M0,M1,M2,M3);
%diag(D)

[V,D] = generalizedEigenSolver(M0,M1,M2,M3);

%diag(D)
V=V;
D=diag(D);
s=D./(D+1);
s=s';
%V=V'
t=V(2,:)./(V(1,:)+V(2,:));

function [V,D] = generalizedEigenSolver(M0,M1,M2,M3)
%M3inv=M3^(-1);
Z=zeros(3,3);
I=eye(3);
C2=[Z I Z;Z Z I;-M0 -M1 -M2];
C1=[I Z Z;Z I Z;Z Z M3];
%detr=det(C1)
%detr=det(C2)
[V,D] = eig(C2,C1,'qz');



function [V,D] = eigenSolver(M0,M1,M2,M3)
M3inv=M3^(-1);
Z=zeros(3,3);
I=eye(3);
C=[Z I Z;Z Z I;-M3inv*M0 -M3inv*M1 -M3inv*M2];
[V,D] = eig(C);

function M=SubMatrixPolynomialFunction(B1E,B2E,n)
Mx=zeros(3,3);
My=zeros(3,3);
Mw=zeros(3,3);


for i=0:2
for j=0:2
k=min(i,j);
for l=0:k
    %l
    m=i+j+1-l;
    if m<=3
    %    [i j l m]
    Mx(i+1,j+1)=Mx(i+1,j+1)+ImplicifCoefficient(B1E,B2E,l,m)*(B1E(2,m+1)-B1E(2,l+1));
    My(i+1,j+1)=My(i+1,j+1)+ImplicifCoefficient(B1E,B2E,l,m)*(B1E(1,l+1)-B1E(1,m+1));
    Mw(i+1,j+1)=Mw(i+1,j+1)+ImplicifCoefficient(B1E,B2E,l,m)*(B1E(2,l+1)*B1E(1,m+1)-B1E(2,m+1)*B1E(1,l+1));
    
    end
end
end
end


%B2E(1,n+1)*Mx+B2E(2,n+1)*My
%B2E(end,n+1)*Mw
M=nchoosek(3,n)*(B2E(1,n+1)*Mx+B2E(2,n+1)*My+B2E(end,n+1)*Mw );

function ic=ImplicifCoefficient(B1E,B2E,l,m)
ic=nchoosek(3,l)*nchoosek(3,m)*B1E(end,l+1)*B1E(end,m+1);

function plotSpline(B,clr)
np=100000;
u=[0:np]/np;
%dim=size(B,1);
R=position(B,u);
%R=repmat(B(:,1),1,np+1).*repmat((1-u).^3,dim,1) + repmat(B(:,2),1,np+1).*repmat(3*(1-u).^2.*u,dim,1) + repmat(B(:,3),1,np+1).*repmat(3*(1-u).*u.^2,dim,1) + repmat(B(:,4),1,np+1).*repmat(u.^3,dim,1);

hold on
%if size(B,1)==2
plot(B(1,:),B(2,:),[clr '--'])
plot(R(1,:),R(2,:),clr,'Linewidth',1)
%elseif size(B,1)==3
%plot3(B(1,:),B(2,:),B(3,:),[clr '--'])
%plot3(R(1,:),R(2,:),R(3,:),clr,'Linewidth',1)
%end
%axis equal
grid on

function R=position(B,u)
np=length(u);
dim=size(B,1);
%size(repmat(B(:,1),1,np))
%size(repmat((1-u).^3,dim,1))
R=repmat(B(:,1),1,np).*repmat((1-u).^3,dim,1) + repmat(B(:,2),1,np).*repmat(3*(1-u).^2.*u,dim,1) + repmat(B(:,3),1,np).*repmat(3*(1-u).*u.^2,dim,1) + repmat(B(:,4),1,np).*repmat(u.^3,dim,1);
%R=B(:,1)*(1-u)^3 + B(:,2)*3*(1-u)^2*u + B(:,3)*3*(1-u)*u^2 + B(:,4)*u^3;
