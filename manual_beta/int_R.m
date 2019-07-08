function int_R
%syms a b c u real

%f=sqrt(a*u^2+b*u+c);

%I=int(f,u)
%Id=simple(subs(I,u,1)-subs(I,u,-1));
%latex(Id)

A=rand(3,1);
B=rand(3,1);
x=rand(3,1);

a=dot(B-A,B-A);
c=dot(A-x,A-x);
b=2*dot(A-x,B-A);

4*a*c-b^2