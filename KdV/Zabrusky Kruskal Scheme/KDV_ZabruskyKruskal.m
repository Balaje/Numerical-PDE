% Zabrusky Kruskal Scheme
% for KdV
% Periodic boundary condition
close all
clear
clc

A = 8;
m = sqrt(A/2);
a = -20; b = 20;
h = 0.1739;
N = fix((b-a)/h);

tf = 4; t0 = 0;
k = 0.002;
M = fix((tf-t0)/k);

x = a:h:b;

u0 = A*sech(m*x).^2;

U1 = zeros(N+1,1);
U2 = U1;

mu = 2*k/h;
tau = k/(h^3);

U1(1) = u0(1)-0.5*mu*(u0(2)+u0(1)+u0(N))*(u0(2)-u0(N))-0.5*tau*(u0(3)-2*u0(2)+2*u0(N)-u0(N-1));
U1(2) = u0(2)-0.5*mu*(u0(3)+u0(2)+u0(1))*(u0(3)-u0(1))-0.5*tau*(u0(4)-2*u0(3)+2*u0(1)-u0(N));
for j=3:N-1
    U1(j) = u0(j)-0.5*mu*(u0(j+1)+u0(j)+u0(j-1))*(u0(j+1)-u0(j-1))-0.5*tau*(u0(j+2)-2*u0(j+1)+2*u0(j-1)-u0(j-2));
end
U1(N) = u0(N)-0.5*mu*(u0(N+1)+u0(N)+u0(N-1))*(u0(N+1)-u0(N-1))-0.5*tau*(u0(2)-2*u0(N+1)+2*u0(N-1)-u0(N-2));
U1(N+1) = u0(N+1)-0.5*mu*(u0(2)+u0(N+1)+u0(N))*(u0(2)-u0(N))-0.5*tau*(u0(3)-2*u0(2)+2*u0(N)-u0(N-1));

for p=2:M
    U2(1) = u0(1)-mu*(U1(2)+U1(1)+U1(N))*(U1(2)-U1(N))-tau*(U1(3)-2*U1(2)+2*U1(N)-U1(N-1));
    U2(2) = u0(2)-mu*(U1(3)+U1(2)+U1(1))*(U1(3)-U1(1))-tau*(U1(4)-2*U1(3)+2*U1(1)-U1(N));
    for j=3:N-1
        U2(j) = u0(j)-mu*(U1(j+1)+U1(j)+U1(j-1))*(U1(j+1)-U1(j-1))-tau*(U1(j+2)-2*U1(j+1)+2*U1(j-1)-U1(j-2));
    end
    U2(N) = u0(N)-mu*(U1(N+1)+U1(N)+U1(N-1))*(U1(N+1)-U1(N-1))-tau*(U1(2)-2*U1(N+1)+2*U1(N-1)-U1(N-2));
    U2(N+1) = u0(N+1)-mu*(U1(2)+U1(N+1)+U1(N))*(U1(2)-U1(N))-tau*(U1(3)-2*U1(2)+2*U1(N)-U1(N-1));
    
    u0 = U1;
    U1 = U2;
    
    plot(x,U1,'b');
    axis([a,b,-A,2*A])
    pause(0.001);
    %hold on
end