% Problem 1 Lab sheet 2 Crank Nicolson Scheme

clc
format long
close 
clear

a = 0; b = 1;
N = 10;

for p=1:5
h(p) = (b-a)/N;
x = a:h(p):b;

t0 = 0; tf = 1;
M = 100;
k = (tf-t0)/M;

mu = k/(2*h(p)^2);

u0 = zeros(N+1,1);
for j=1:N+1
    u0(j) = cos(pi*x(j));
end

A = zeros(N+1);
B = zeros(N+1);
F = zeros(N+1,1);

A(1,1) = 1+2*mu;
B(1,1) = 1-2*mu;
for j=2:N+1
    A(j,j) = 1+2*mu;
    B(j,j) = 1-2*mu;
    A(j-1,j) = -mu;
    A(j,j-1) = -mu;
    B(j-1,j) = mu;
    B(j,j-1) = mu;
end
A(1,2) = -2*mu;
A(N+1,N) = -2*mu;
B(1,2) = 2*mu;
B(N+1,N) = 2*mu;

u = zeros(N+1,1);
t = t0;
for i=1:M
    for j=1:N+1
        F(j) = 0.5*k*(f1(x(j),t)+f1(x(j),t+k));
    end
    u = A\(B*u0+F);
    
    u0 = u;
    t = t+k;
end

exact = zeros(N+1,1);
for j=1:N+1
    exact(j) = exp(-tf)*cos(pi*x(j));
end

error(p) = max(abs(u-exact));
N = N*2;
end

for j=1:p-1
    order(j) = log(error(j)/error(j+1))/log(2);
end
plot(x,u,x,exact);
