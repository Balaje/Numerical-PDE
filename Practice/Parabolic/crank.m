% Program to implement Crank Nicolson Method for
% Problem 2 lab-sheet 2
clear 
close all

%format long

N = 10;
M = 1000;

for q=1:5

a = 0; b = 1;
h(q) = (b-a)/N;

t0 = 0; tf = 1;
k = (tf-t0)/M;

mu = k/(2*h(q)^2);

x = a:h(q):b;

u0 = zeros(N+1,1);
for j=1:N+1
    u0(j) = x(j)*(1-x(j));
end

A = zeros(N+1);
B = zeros(N+1);
u = zeros(N+1,1);

A(1,1) = 1+2*mu-2*mu*h(q);
for j=2:N+1
    A(j,j) = 1+2*mu;
    A(j,j-1) = -mu;
    A(j-1,j) = -mu;
end
A(1,2) = -2*mu;
A(N+1,N) = -2*mu;

B(1,1) = 1-2*mu+2*mu*h(q);
for j=2:N+1
    B(j,j) = 1-2*mu;
    B(j,j-1) = mu;
    B(j-1,j) = mu;
end
B(1,2) = 2*mu;
B(N+1,N) = 2*mu;

F = zeros(N+1,1);
t = t0;

for p=1:M-1
    F(1) = k*0.5*(f(x(1),t)+f(x(1),t+k)) - 2*mu*h(q)*exp(-t) - 2*mu*h(q)*exp(-(t+k));
    for j=2:N
        F(j) = k*0.5*(f(x(j),t)+f(x(j),t+k));
    end
    F(N+1) = k*0.5*(f(x(N+1),t)+f(x(N+1),t+k)) - 2*mu*h(q)*exp(-t)- 2*mu*h(q)*exp(-(t+k));
    
    u = A\(B*u0+F);
    
    u0 = u;
    t = t+k;
end
exact = zeros(N+1,1);
for j=1:N+1
    exact(j) = exp(-tf)*(x(j)*(1-x(j)));
end

error(q) = max(abs(u-exact));

N = N*2;
end

for j=1:q-1
    order(j) = log(error(j)/error(j+1))/log(2);
end
figure(1)
plot(x,u,x,exact);
figure(2)
plot(log(h),log(error));