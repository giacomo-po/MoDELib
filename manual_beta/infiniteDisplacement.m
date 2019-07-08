%function infiniteDisplacement
syms Ax Ay Az Bx By Bz qx qy qz u L x y z real

A=[Ax Ay Az]';        % initial point of the source segment
B=[Bx By Bz]';        % final point of the source segment
t=(B-A)/L;            % the tangent vector
S=A*(1-u)+B*u;        % the source point moving along A->B as a function of u in [0,1]
X=[x y z]';           % the field point
Rv=X-S;               % the R vector (source to field)
R=sqrt(dot(Rv,Rv));   % the magnitude of the R vector
q=[qx qy qz]';        % The q vector


%I=int(dot(t,cross(q,Rv))/(R*(R-dot(q,Rv))),u,0,1);
I=int(Rv/dot(Rv,Rv));
disp('Done integration')
Is=simple(I);
Is


