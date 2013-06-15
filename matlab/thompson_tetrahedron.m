function [N,V,R]=thompson_tetrahedron(L)
%clc
%clear
%close all

alpha=1;
beta=2;
gamma=3;
delta=4;


N(:,delta) = [ 1.0, 1.0, 1.0]';		% delta
N(:,beta)  = [ 1.0,-1.0,-1.0]';		% beta
N(:,alpha) = [-1.0, 1.0,-1.0]';		% alpha
N(:,gamma) = [-1.0,-1.0, 1.0]';		% gamma

%N=N+repmat([0 0 1]',1,4);

planenames={'alpha', 'beta', 'gamma', 'delta'};


N=N/sqrt(3)*L;


A=alpha;
B=beta;
C=gamma;
D=delta;

V(:,A) = -3*N(:,delta);		% delta
V(:,B)  = -3*N(:,beta);		% beta
V(:,C) = -3*N(:,alpha);		% alpha
V(:,D) = -3*N(:,gamma);		% gamma


figure(gcf)
%clf
hold on

%plot3(V(1,:),V(2,:),V(3,:),'ok','Linewidth',2)
%plot3(0,0,0,'og','Linewidth',2)
clr='cryo';
for k=1:4
  %  quiver3(0,0,0,N(1,k),N(2,k),N(3,k),0,'g','Linewidth',2)
    text(N(1,k),N(2,k),N(3,k),texlabel(planenames{k}),'Fontsize',14,'Color','r')
   % quiver3(0,0,0,V(1,k),V(2,k),V(3,k),0,'m','Linewidth',2)
    ka=rem(k+[0:2],4)+1;
    patch(V(1,ka),V(2,ka),V(3,ka),clr(k),'Facealpha',0.3)
end
axis equal
grid on
xlabel('x')
ylabel('y')
zlabel('z')

for k=1:4
    %k
R(:,3,k)=N(:,k)/norm(N(:,k));
k1=rem(k,4)+1
v1=cross(R(:,3,k),N(:,k1)/norm(N(:,k1)));
R(:,1,k)=v1/norm(v1);
R(:,2,k)=cross(R(:,3,k),R(:,1,k));
quiver3(N(1,k),N(2,k),N(3,k),R(1,1,k)*L,R(2,1,k)*L,R(3,1,k)*L,'k','Linewidth',1)
quiver3(N(1,k),N(2,k),N(3,k),R(1,2,k)*L,R(2,2,k)*L,R(3,2,k)*L,'r','Linewidth',1)
quiver3(N(1,k),N(2,k),N(3,k),R(1,3,k)*L,R(2,3,k)*L,R(3,3,k)*L,'g','Linewidth',1)

end


