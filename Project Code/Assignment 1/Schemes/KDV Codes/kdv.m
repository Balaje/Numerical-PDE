clear
close all
clc
format long

a = -20; b = 20;
t0 = 0; tf = 10;

N = 600;

h = (b-a)/N;
tau = 0.2;

k = tau*h^3;

mu = 3*k/h;

x = a+h:h:b-h;
u0 = zeros(N-1,1);
un = u0;

u0 = 8*(sech(2*(x+8))).^2 + 4.5*(sech(3*(x+5))).^2;
%plot(x,u0)
%u0 = cos(pi*x);
M = fix((tf-t0)/k);

for i=1:M
    un(1) = u0(1) - mu*(f(u0(2)) - f(0)) - tau*(u0(3)-3*u0(2)+3*u0(1));
    for j=2:N-3
        un(j) = u0(j) - mu*(f(u0(j+1)) - f(u0(j-1))) - tau*(u0(j+2)-3*u0(j+1)+3*u0(j)-u0(j-1));
    end
    un(N-2) = u0(N-2) - mu*(f(u0(N-1)) - f(u0(N-3))) - tau*(-3*u0(N-1)+3*u0(N-2)-u0(N-3));
    un(N-1) = u0(N-1) - mu*(f(0) - f(u0(N-2))) - tau*(3*u0(N-1)-u0(N-2));
    %str = ['Time: ', num2str(i*k)];
    
    
    plot(x,un);
    axis([a,b,-1,8]);
    grid on
    %hold on;
    %title(str)
    pause(0.01)
    u0 = un;
end

%plot(x,un)