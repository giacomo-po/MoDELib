function plotLattice(A,R,marker,linewidth)

hold on
i=0;
j=0;
k=0;

while norm(A*[i j k]')<R
    while norm(A*[i j k]')<R
        while norm(A*[i j k]')<R
            P0=A*[i j k]';
            plot3(P0(1),P0(2),P0(3),marker,'Linewidth',linewidth)
            k=k+1;
        end
        k=0;
        while norm(A*[i j k]')<R
            P0=A*[i j k]';
            plot3(P0(1),P0(2),P0(3),marker,'Linewidth',linewidth)
            k=k-1;
        end
        k=0;
        j=j+1;
    end
    j=0;
    while norm(A*[i j k]')<R
        while norm(A*[i j k]')<R
            P0=A*[i j k]';
            plot3(P0(1),P0(2),P0(3),marker,'Linewidth',linewidth)
            k=k+1;
        end
        k=0;
        while norm(A*[i j k]')<R
            P0=A*[i j k]';
            plot3(P0(1),P0(2),P0(3),marker,'Linewidth',linewidth)
            k=k-1;
        end
        k=0;
        j=j-1;
    end
    j=0;
    i=i+1;
end
i=0;
while norm(A*[i j k]')<R
    while norm(A*[i j k]')<R
        while norm(A*[i j k]')<R
            P0=A*[i j k]';
            plot3(P0(1),P0(2),P0(3),marker,'Linewidth',linewidth)
            k=k+1;
        end
        k=0;
        while norm(A*[i j k]')<R
            P0=A*[i j k]';
            plot3(P0(1),P0(2),P0(3),marker,'Linewidth',linewidth)
            k=k-1;
        end
        k=0;
        j=j+1;
    end
    j=0;
    while norm(A*[i j k]')<R
        while norm(A*[i j k]')<R
            P0=A*[i j k]';
            plot3(P0(1),P0(2),P0(3),marker,'Linewidth',linewidth)
            k=k+1;
        end
        k=0;
        while norm(A*[i j k]')<R
            P0=A*[i j k]';
            plot3(P0(1),P0(2),P0(3),marker,'Linewidth',linewidth)
            k=k-1;
        end
        k=0;
        j=j-1;
    end
    j=0;
    i=i-1;
end


end

