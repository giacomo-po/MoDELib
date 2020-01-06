clear all
close all
clc

syms A1 A2 A3 B1 B2 B3 X1 X2 X3 s1 s2 s3 u real

A=[A1 A2 A3]';
B=[B1 B2 B3]';
X=[X1 X2 X3]';

s=[s1 s2 s3]';
X1=A+u*(B-A)

R=X-X1;

f=simplify(R/norm(R)/(norm(R)+dot(s,R)))
I=int(f,u,0,1)
