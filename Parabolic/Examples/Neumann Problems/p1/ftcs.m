clear all;
close all;

format long;

b = 1; a = 0;
t0 = 0; tf = 1;
N = 15;

h = (b-a)/N;

mu = 0.2;

k = mu*h^2;

M = fix((tf-t0)/k);

x = a:h:b;
u0 = zeros(N+1,1);
u_new = u0;

for j=1:N+1
    u0(j) = cos(pi*x(j));
end

t=0;
for i=1:M-1
    u_new(1) = u0(1) + mu*(u0(2) - 2*u0(1) + u0(2)) + k*f(t,x(1));
    for j=2:N
        u_new(j) = u0(j) + mu*(u0(j-1) - 2*u0(j) + u0(j+1))+ k*f(t,x(j));
    end
    u_new(N+1) = u0(N+1) + mu*(u0(N) - 2*u0(N+1) + u0(N)) + k*f(t,x(N+1));
    u0 = u_new;
    t=t+k;
end

exact = zeros(N+1,1);
for j=1:N+1
    exact(j) = exp(-t)*cos(pi*x(j));
end

plot(x,u_new,'b--',x,exact,'g-+');