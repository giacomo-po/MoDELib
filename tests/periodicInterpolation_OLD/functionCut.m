function f=functionCut(data,d,L,np)
d=d/norm(d);
dL=L/(np-1);
f=zeros(np,1);

for k=1:np
x=(k-1)*dL*d;
f(k)=0;
for i=1:size(data,1)
    kv=data(i,[1 2]);
    S=data(i,3);
    C=data(i,4);
    f(k)=f(k)+S*sin(dot(kv,x))+C*cos(dot(kv,x));
end
end

end