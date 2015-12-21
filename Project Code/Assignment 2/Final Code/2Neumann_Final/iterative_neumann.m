clear all;
close all;
clc

a=0;b=2;
c=0;d=1;
N=100;
hx=(b-a)/N;

hx1=hx*hx;
hy=hx;
hy1=hx1;

M = fix((d-c)/hx);

x=a:hx:b;
y=c:hy:d;
tic
u1=zeros(N+1,M+1);
u2=zeros(N+1,M+1);
usol=zeros(N+1,M+1);

for i=1:N+1
    u2(i,1)=1;
    u2(i,M+1)=exp(x(i));
end

for j=1:M+1
    u2(1,j)=1;
    u2(N+1,j)=exp(2*y(j));
end

tol = 1E-6;

while max(max(abs(u2-u1)))>tol
    
    u1=u2;
    
    for i=2:N
        for j=2:M
            u2(i,j) = ((hx1*hy1)/(2*hx1+2*hy1))*((u1(i+1,j)+u1(i-1,j))/hx1 + (u1(i,j+1)+u1(i,j-1))/hy1 - f(x(i),y(j)));
        end
    end
end

figure(1)
surf(y,x,u2)
title('The approximate solution');

for i=1:N+1
    for j=1:M+1
        usol(i,j)=exp(x(i)*y(j));
    end
end

figure(2);
surf(y,x,usol);
title('The exact solution')
toc