% Two way wave equation solution
% u_tt = u_xx
% I.C : u(x,0) = 1/8(sin(pi*x)) u_t(x,0) = 0
% B.C : u(0,t) = u(1,t) = 0
clc
clear
close all
format long

a = 0; b = 1;
N = 100;

for p=1:7
h(p) = (b-a)/N;

t0 = 0; tf = 1;
k = 0.5*h(p);
M = fix((tf-t0)/k);
mu = (k/h(p))^2;

x = a+h(p):h(p):b-h(p);
u0 = zeros(N-1,1);
for j=1:N-1
    u0(j) = 1/8*sin(pi*x(j));
end

u1 = zeros(N-1,1);
u = zeros(N-1,1);
% n = 1
u1(1) = 0.5*( 2*u0(1) + mu*(u0(2)-2*u0(1)) );
for j=2:N-2
    u1(j) = 0.5*( 2*u0(j) + mu*(u0(j-1)-2*u0(j)+u0(j+1)) );
end
u1(N-1) = 0.5*( 2*u0(N-1) + mu*(u0(N-2)-2*u0(N-1)) );

% n = 2:M
t=t0;
for i=2:M
    u(1) = 2*u1(1) - u0(1) + mu*(u1(2)-2*u1(1));
    for j=2:N-2
        u(j) = 2*u1(j) - u0(j) + mu*(u1(j-1)-2*u1(j)+u1(j+1));
    end
    u(N-1) = 2*u1(N-1) - u0(N-1) + mu*(u1(N-2)-2*u1(N-1));
    
    u0 = u1;
    u1 = u;
    t = t+k;
end

ue = zeros(N-1,1);
for j=1:N-1
    ue(j) = 1/8*sin(pi*x(j))*cos(pi*t);
end

error(p) = max(abs(u-ue));

N = 2*N;
end

for j=1:p-1
    order(j) = log(error(j)/error(j+1))/log(2);
end
plot(x,u,x,ue)