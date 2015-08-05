clear all
close all
clc

syms x y z f real

syms f000 f000x f000y f000z real
syms f100 f100x f100y f100z real
syms f010 f010x f010y f010z real
syms f001 f001x f001y f001z real

dof=[f000 f000x f000y f000z...
     f100 f100x f100y f100z...
     f010 f010x f010y f010z...
     f001 f001x f001y f001z]';

%dof=[dof;]

     f000+[f000x f000y f000z]*([1 1 0]/3-[0 0 0])'...
    +f100+[f100x f100y f100z]*([1 1 0]/3-[1 0 0])'...
    +f010+[f010x f010y f010z]*([1 1 0]/3-[0 1 0])'

return

p=3;
f=0;
Na=0;
for i=0:p
    for j=0:p-i
        for k=0:p-i-j
        eval(['syms a' num2str(i) num2str(j) num2str(k) ' real'])
        eval(['f=f+a' num2str(i) num2str(j) num2str(k) '*x^i*y^j*z^k;'])
        Na=Na+1;
        end
    end
end


%%
df=[diff(f,x);diff(f,y);diff(f,z)];
fdf=[f;df];
lhs=[subs(fdf,{x,y,z},{0,0,0});
    subs(fdf,{x,y,z},{1,0,0});
    subs(fdf,{x,y,z},{0,1,0});
    subs(fdf,{x,y,z},{0,0,1})]
Ne=size(lhs,1);

%C=zeros(Na,Na);
%A=zeros(Na,1);

f011=subs(f,{x,y,z},{0,1/3,1/3});
f101=subs(f,{x,y,z},{1/3,0,1/3});
f110=subs(f,{x,y,z},{1/3,1/3,0});
f111=subs(f,{x,y,z},{1/3,1/3,1/3});

n=1;
for i=0:p
    for j=0:p-i
        for k=0:p-i-j
        eval(['C(:,n)=diff(lhs,a' num2str(i) num2str(j) num2str(k) ');'])
        eval(['C(Ne+1,n)=diff(f011,a' num2str(i) num2str(j) num2str(k) ');']) % value at 011/3 (face center)
        eval(['C(Ne+2,n)=diff(f101,a' num2str(i) num2str(j) num2str(k) ');']) % value at 101/3 (face center)
        eval(['C(Ne+3,n)=diff(f110,a' num2str(i) num2str(j) num2str(k) ');']) % value at 110/3 (face center)
        eval(['C(Ne+4,n)=diff(f111,a' num2str(i) num2str(j) num2str(k) ');']) % value at 111/3 (face center)
        n=n+1;
        end
    end
end

%%
%A=inv(C)*dof
