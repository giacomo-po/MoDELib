close all
clear all
clc

s0=rand(6,1);
e0=rand(6,1);
C=rand(6,6);
C=C'*C;
eP=rand(6,1);

alpha=[0 0 1e10 0 1e5 0];
%alpha=[1 1 1 1 1 1]*1e10;
%alpha=[0 0 0 0 0 0];

strategy=2;

switch strategy
    case 1 % OLD WRONG STRATEGY
        A=diag(alpha);
        A1=diag(1./(1+alpha));
        A2=diag(alpha./(1+alpha));
        s=A1*s0+A2*C*(e0-eP);
    case 2 % NEW WORKING STRATEGY
        A=diag(alpha)*diag(diag(C));
        s=inv(eye(6)+A*inv(C))*(s0+A*(e0-eP));
end
e=inv(C)*s+eP;

[e e0 s s0]
