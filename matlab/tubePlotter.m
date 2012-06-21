function TubePlotter(H,alpha,marker,lw)
%alpha=0.5;
g=parametricLength(H(1,:),H(3,:),alpha);
%g = norm(-H(3,:))^alpha;
np=100;
u=[0:(np-1)]'/(np-1);

UPOW=[u.^0 u u.^2 u.^3];
UPOWu=[u.^0*0 u.^0 2*u 3*u.^2];
UPOWuu=[u.^0*0 u*0 u.^2*0+2 6*u];

SFC=[1    0  0  0;
    0    g  0  0;
    -3 -2*g  3 -g;
    2    g -2  g];
% UPOW*SFC

R=UPOW*SFC*H;
T=UPOWu*SFC*H;
%size(R)
plot3(R(:,1),R(:,2),R(:,3),marker,'Linewidth',lw)
plot3(R(1,1),R(1,2),R(1,3),'.r')
plot3(R(end,1),R(end,2),R(end,3),'.r')
quiver3(R(end,1),R(end,2),R(end,3),T(end,1),T(end,2),T(end,3),'.m')
quiver3(R(1,1),R(1,2),R(1,3),T(1,1),T(1,2),T(1,3),'.m')

return

T=UPOWu*SFC*H;
%K=UPOWuu*SFC*H/g^2;
%size(R)

%dim=size(B,1);
%R=position(B,u);
%R=repmat(B(:,1),1,np+1).*repmat((1-u).^3,dim,1) + repmat(B(:,2),1,np+1).*repmat(3*(1-u).^2.*u,dim,1) + repmat(B(:,3),1,np+1).*repmat(3*(1-u).*u.^2,dim,1) + repmat(B(:,4),1,np+1).*repmat(u.^3,dim,1);

%hold on
%if size(B,1)==2
%plot(B(1,:),B(2,:),[clr '--'])
%plot(R(1,:),R(2,:),clr,'Linewidth',2)
%elseif size(B,1)==3
npol=10;
radius=0.08;
FaceAlpha=0.3;
EdgeAlpha=0.1;
fr=1.1;
[sx,sy,sz]=sphere(npol);
[cx,cy,cz]=cylinder(radius*[1:np],npol);



figure(1)
hold on
%% Plot spheres at the nodes

surf(sx*fr*radius+H(1,1),sy*fr*radius+H(1,2),sz*fr*radius+H(1,3),sz*0,'FaceColor','k','EdgeColor','k','FaceAlpha',FaceAlpha,'EdgeAlpha',EdgeAlpha);
surf(sx*fr*radius+H(3,1),sy*fr*radius+H(3,2),sz*fr*radius+H(3,3),sz*0,'FaceColor','k','EdgeColor','k','FaceAlpha',FaceAlpha,'EdgeAlpha',EdgeAlpha);

%% 
%plot3(R(:,1),R(:,2),R(:,3),'b')

%plot3(H(1,1),H(1,2),H(1,3),'xg')
%plot3(H(3,1),H(3,2),H(3,3),'xg')

%np=size(R,1);


%cxr=reshape(cx,prod(size(cx)),1)';
%cyr=reshape(cy,prod(size(cy)),1)';
%czr=reshape(cz,prod(size(cz)),1)';
temp=rand(3,1);
for k=1:np
Q3=T(k,:)'/norm(T(k,:));

Q1=cross(Q3,temp);
Q1=Q1/norm(Q1);
Q2=cross(Q3,Q1);
Q=[Q1 Q2 Q3];
%det(Q)

CR=Q*[cx(1,:);cy(1,:);cz(1,:)];

CRx(k,:)=CR(1,:)+R(k,1);
CRy(k,:)=CR(2,:)+R(k,2);
CRz(k,:)=CR(3,:)+R(k,3);
%CRy=reshape(CR(2,:),2,size(CR,2)/2);
%CRz=reshape(CR(3,:),2,size(CR,2)/2);

end
h= surf(CRx,CRy,CRz,CRz*0,'FaceColor','b','EdgeColor','b','FaceAlpha',FaceAlpha,'EdgeAlpha',EdgeAlpha);

%set(h,'FaceLighting','flat','FaceColor','interp','AmbientStrength',0.5)
%light('Position',[1 0 0],'Style','infinite')

return
%% shade area inside
x=[R(:,1); R(end,1); R(1,1)];
y=[R(:,2); 0; 0];
z=[R(:,3); 0; 0];
hold on
plt1=patch(x,y,z,0,'FaceColor','k','FaceAlpha',0.2,'EdgeAlpha',0.0);


