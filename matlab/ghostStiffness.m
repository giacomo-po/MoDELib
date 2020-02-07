clear all
close all
clc

syms f L real
K0=[L/3 L/6;
    L/6 L/3];

K1=K0*f;        % stiffness matrix of segment of length fL
K2=K0*(1-f);    % stiffness matrix of segment of length (1-f)L  

K=[K1(1,1) K1(1,2) 0;K1(2,1) K1(2,2)+K2(1,1) K2(1,2);0 K2(2,1) K2(2,2)] % stiffness matrix of the two segments
M=[1 0;1-f f;0 1] % impose that midpoint velocity is interpolated linearly

Kr=simplify(M'*K*M) % reduced stiffness withou midpont node
