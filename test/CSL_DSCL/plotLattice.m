function plotLattice(A,L,marker,linewidth)

hold on
invA=inv(A);

P=[-1  1  1  1 -1  1 1  1;
    -1 -1  1 -1 -1 -1 1 -1;
    -1 -1 -1  0  1  1 1  1]*L;
N=round(invA*P);

% iRange=[min(N(1,:))-1:max(N(1,:))+1];
% jRange=[min(N(2,:))-1:max(N(2,:))+1];
% kRange=[min(N(3,:))-1:max(N(3,:))+1];

iRange=[min(min(N))-1:max(max(N))+1];
jRange=iRange;
kRange=iRange;

tol=1e-5;

X=[];
for i=iRange
    for j=jRange
        for k=kRange
            P0=A*[i j k]';
            if (max(P0)<L+tol && min(P0)>-L-tol)
                X=[X P0];
            end
            
        end
    end
end
plot3(X(1,:),X(2,:),X(3,:),marker,'Linewidth',linewidth)

end

