function C2G=integerc2G(R)
for r=1:size(R,1)
X=R(r,:);
nMax=max(abs(X));
[N,D] = rat(X/nMax);
d = lcm(D(1),lcm(D(2),D(3)));
%B=inv(A');
C2G(r,:)=d*N./D;
%r=B*R;
end