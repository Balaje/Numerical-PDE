clear
close all
clc
format long

a = 0; b = 1;
t0 = 0; tf = 0.3;

N = 320;

h = (b-a)/N;
tau = 0.2;

k = tau*h^3;

mu = 6*k/h;

x = a:h:b;
u0 = zeros(N+1,1);
un = u0;

%u0 = 12.5*(sech(2.5*(x+2))).^2 + 8*(sech(2*(x+1))).^2;
u0 = cos(pi*x);
M = fix((tf-t0)/k);

for i=1:M
    un(N+1) = u0(1);
    un(1) = u0(1) - mu*(f(u0(2)) - f(u0(1))) - tau*(u0(4)-3*u0(3)+3*u0(2)-u0(1));
    for j=2:N-3
        un(j) = u0(j) - mu*(f()-f())
    end
    
end

plot(x,un)
