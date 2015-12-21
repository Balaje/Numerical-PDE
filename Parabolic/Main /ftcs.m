clear all;
close all;

clc

format long;

N = 10; % Number of subintervals
a = 0; b = 1; 
h = (b-a)/N;
alpha = 1;

% Defining mu directly
mu = 0.2;

% Using this define time step.
del_t = mu*h^2/alpha;

% Define initial and final times.
t0 = 0; tf = 1;
M = fix((tf-t0)/del_t);

x = zeros(N-1,1);
u0 = zeros(N-1,1);
u_new = u0;

t = t0:del_t:tf;

for i=1:N-1
    x(i) = i*h;
    u0(i) = sin(pi*x(i));
end

for i=2:M
    u_new(1) = (1-2*mu)*u0(1) + mu*u0(2) + del_t*(f(t(i),x(1)));
    u_new(N-1) = mu*u0(N-2) + (1-2*mu)*u0(N-1) + del_t*(f(t(i),x(N-1)));
    for j=2:N-2
        u_new(j) = mu*u0(j-1)+ (1-2*mu)*u0(j) + mu*u0(j+1) + del_t*f(t(i),x(j));
    end
    u0 = u_new;
end

% Compare with the exact solution
exact = zeros(N-1,1);
for j=1:N-1
    exact(j) = exp(-t(i))*sin(pi*x(j));
end

plot(x,u_new,'b--',x,exact,'r-+');