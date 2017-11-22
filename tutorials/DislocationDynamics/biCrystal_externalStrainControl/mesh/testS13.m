

X=[4 -3 -1];
X=X/norm(X);

Y=[1 1 1];
Y=Y/norm(Y);

Z=[-2 -5 7];
Z=Z/norm(Z);

R=[X;Y;Z]

[V,D]=eig(R)