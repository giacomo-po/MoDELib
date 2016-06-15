clear all
close all
clc

N=load('N/N_6.txt');
K=load('K.txt');


load MA.mat;
sol=load('sol.txt');

% Find permutation matrix used by COMSOL
for(k=1:size(sol,1))
nd=sum(((N(:,[2:4])-repmat(sol(k,[1:3]),size(sol,1),1)).^2)');
P(k,1)=find(nd==min(nd));
end
if (length(unique(P))~=length(P))
error('permutation not found')
end

posNorm=norm(N(P,[2 3 4])-sol(:,[1 2 3]))

%K=K(P,P);
[V,D]=eig(K);
d=diag(D)
condK=cond(K)


Kc=full(MA.K);
[Vc,Dc]=eig(Kc);
dc=diag(Dc);
condKc=cond(Kc)

eigenNorm=norm(d-dc)

for(k=1:size(sol,1))
nd=abs(V(:,1)-repmat(Vc(k,1),size(Vc,1),1));
%min(nd)
P1(k,1)=find(nd==min(nd));
end
if (length(unique(P1))~=length(P1))
error('permutation not found')
end

K=K(P1,P1); % apply permutation to K
N=N(P,:); % apply permutation to N

normKdiff=norm(K-Kc)/norm(Kc)



%return

%x=linsolve(Kc,MA.L);
%plot(x/max(abs(x)),'r')
b=zeros(size(K,1),1);
center=sum((sol(:,[1 2 3]).^2)');
%return
center=find(center==min(center))
%return
b(center)=1;
x1=K\b;

figure(1)
hold on
plot(sol(:,4)/max(abs(sol(:,4))),'k--')
plot(x1/max(abs(x1)),'g-.')

%comsolNorm=norm(Kc*x-MA.L)

figure
%scatter3(V(:,2),V(:,3),V(:,4),10,MA.L)
circSizes=sol(:,4)-min(sol(:,4));
circSizes=circSizes/max(circSizes);
scatter3(sol(:,1),sol(:,2),sol(:,3),100*(circSizes+0.1),sol(:,4))

figure
circSizes=x1-min(x1);
circSizes=circSizes/max(circSizes);
scatter3(N(:,2),N(:,3),N(:,4),100*(circSizes+0.1),x1)