function r=reciprocalLatticeDirection(n,A)
X=A'*n;
nMax=max(abs(X));
[N,D] = rat(X/nMax);
d = lcm(D(1),lcm(D(2),D(3)));
B=inv(A');
R=d*N./D;
r=B*R;