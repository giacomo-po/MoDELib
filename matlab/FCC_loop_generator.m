%function FCC_loop_generator
clc
clear
close all


figure(1)
clf
hold on

r=25

L=0.75*r
[N,V,R]=thompson_tetrahedron(L);

np=30;
theta=[0:np-1]/np*2*pi+pi/6;
%theta=[0:np-1]/np*2*pi+rand(1);
XL=[r*cos(theta);r*sin(theta);zeros(1,np)];

%repmat(N(:,1),1,np)
%k=1
%XL=XL+repmat(N(:,1),1,np)


file_1 = fopen('DDinput.txt','w')
%file_2 = fopen('output2.txt','w')

%nodeformat='node = %4.10f %4.10f %4.10f; \n'
nodeformat='node = %6.15d %6.15d %6.15d; \n'
linkformat='link = %i %i; \n'
constrType='constraintType = %i; \n'
constrformat='constraintNormal = %i %i %i; \n'
flowformat='flow = %i %i %i; \n'
%fclose(file_1)
%fclose(file_2)

%for k=4%:4
%R(:,3)=N(:,k)/norm(N(:,k));
%v1=cross(R(:,3),rand(3,1));
%R(:,1)=v1/norm(v1);
%R(:,2)=cross(R(:,3),R(:,1));

loop=1
k=4 % the 111 plane
s1= [1 0 0];
f1= [1 0 0];
%XG=R(:,:,k)*XL+repmat(N(:,k)/norm(N(:,k))*100,1,np);
%hold on
LL=100;
for k=1:np
XG(:,k)=LL*[0 .8 .8]'+(k-1)*s1'*LL/(np-1);
end
plot3(XG(1,:),XG(2,:),XG(3,:),'k')
text(XG(1,:),XG(2,:),XG(3,:),num2str([0:np-1]'));
%axis equal

%return

%disp(['// Loop ' ])

for n=1:np
fprintf(file_1,nodeformat, XG(:,n)')
fprintf(file_1,constrType, 0 )
fprintf(file_1,constrformat, [0 0 0] )
end

for n=1:np-1
    fprintf(file_1,linkformat, (loop-1)*np+[n-1 n] )
    fprintf(file_1,flowformat, f1 )
end
%fprintf(file_1,linkformat, (loop-1)*np+[np-1 0] )
%fprintf(file_1,flowformat, f1)
xlabel('x')
ylabel('y')
zlabel('z')

return









