%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%          Iterative method for Laplace Equation      %%%%
%%%     u_{xx} + u_{yy} = 0, 0 < x < 1, 0 < y < 1        %%%%
%%%  u(x, 0) = 0, u(0, y) = 0, u(x, 1) = y, u(1, y) = y  %%%%  
%%%           Exact soln: u(x, y) = x*y                  %%%%
%%%        Course: MATH F422, Dr. P. Dhanumjaya          %%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; 
clear all;
N = 20;
tol = 1e-6;
a = 0;
b = 1;
c = 0;
d = 1;
h = (b-a)/N;
h1 = h*h;
x = a:h:b;
y = c:h:d;
u2 = zeros(N+1,N+1);
usol=zeros(N+1,N+1);
omega = 1.5;

for i = 1:N+1
    u2(1,i) = 0; 
    u2(N+1,i) = y(i); 
    u2(i,1) = 0; 
    u2(i,N+1) = x(i); 
end

err = 1000;
k = 0;
u = u2;
while err > tol
    for i = 2:N
        for j = 2:N
            u(i,j) = (1-omega)*u2(i,j) + omega*(u(i-1,j) + u(i+1,j) ...
                + u(i,j-1) + u(i,j+1))/4;
        end
       
    end
    err = max(max(abs(u-u2)));
    u2 = u;
    k = k+1;
end

figure(1);
mesh(x,y,u2);
title('The approximate solution')

for i=1:N+1
    for j=1:N+1
        usol(i,j)=x(i)*y(j);
    end
end

figure(2);
mesh(x,y,usol);
title('The exact solution')

Error = max(max(abs(u-usol)))
