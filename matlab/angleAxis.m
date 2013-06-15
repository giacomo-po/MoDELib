function R=angleAxis(u,theta)

if size(u)==[3 1]
    u=u/norm(u);
R=eye(3)*cos(theta)+sin(theta)*[0 -u(3) u(2);u(3) 0 -u(1);-u(2) u(1) 0]+(1-cos(theta))*u*u';
elseif size(u)==[3 1]
R=angleAxis(u');
else
    error('SIZE U MUST BE [3 1] or [1 3]')
end
    
    

