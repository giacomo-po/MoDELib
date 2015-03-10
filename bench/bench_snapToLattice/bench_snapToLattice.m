clear all
close all
clc

A=[1 0 1;
   1 1 0;
   0 1 1]/2;

invA=inv(A)

N =[-1    1   -1    1;
     1   -1   -1    1;
    -1   -1    1    1];

n=invA*N

D=cross(N(:,1),N(:,2));
d=invA*D
gd=gcd(d(1),gcd(d(2),d(3)));
d=d/gd
D=A*d
