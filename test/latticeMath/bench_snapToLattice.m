clear all
close all
clc

crystal=2

switch crystal
    case 1 % BCC
        A=[  1 -1  1;
             1  1 -1;
            -1  1  1]/sqrt(3); % lattice vectors in units of b
        
        N =[ 1 -1   1 -1  0  0;
             1  1   0  0  1 -1;
             0  0   1  1  1  1];
        
    case 2 % FCC
        A=[1 0 1;
           1 1 0;
           0 1 1]*sqrt(2)/2; % lattice vectors in units of b
        
        N =[-1    1   -1    1;
             1   -1   -1    1;
            -1   -1    1    1]*sqrt(2);
end

%% compute direction
invA=inv(A);
cofA=invA*det(A);

n=invA*N

N1=N(:,2);
N2=N(:,3);
n1=invA*N1;
n2=invA*N2;

D=cross(N1,N2) % common direction in global cordinates

D1=cofA'*cross(n1,n2)

return


d=invA*D % common direction in lattice cordinates
gd=gcd(d(1),gcd(d(2),d(3)));
d=d/gd
D=A*d
