function BZ=H2BZplanar(H,alpha)
g=parametricLength(H(1,:),H(3,:),alpha);
BZ=[H(1,:);
    H(1,:)+H(2,:)/3*g;
    H(3,:)-H(4,:)/3*g;
    H(3,:)];
BZ=BZ(:,[1 2])';

Hp=H;
Hp(:,3)=0;
%tubeplotter(Hp,alpha,'k-')