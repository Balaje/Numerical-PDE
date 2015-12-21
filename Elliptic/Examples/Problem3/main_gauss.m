% Problem 5 Elliptic PDE Labsheet

a = 0;
b = 1;
c = 0;
d = 1;

M = 40; % x
N = 80; % y

hx = (b-a)/M;
hx1 = hx*hx;
hy = (d-c)/N;
hy1 = hy*hy;

x = a:hx:b;
y = c:hy:d;

u2 = zeros(M+1,N+1);
for j=1:N+1
    u2(1,j) = 0.5*cos(2*pi*y(j));
    u2(M+1,j) = 0.5*cos(2*pi*y(j));
end
tol = 10^-6;
error = 1000;
k=0;
u = u2;
%u(i,j) = 0.5*(hx1+hy1)*( (u(i-1,j) + u(i+1,j))/hx1 + 2*(u(i,j))/hy1 - f(x(i),y(j)) );
while error > tol
    for i=2:M
        u(i,1) = (2/hx1+2/hy1)^-1*( (u(i-1,1) + u(i+1,1))/hx1 + (2*u(i,2))/hy1 - f(x(i),y(2)) );
        for j=2:N
           u(i,j) = (2/hx1+2/hy1)^-1*( ( u(i-1,j) + u(i+1,j) )/hx1 + ( u(i,j+1) + u(i,j-1) )/hy1 - f(x(i),y(j)) ); 
        end
        u(i,N+1) = (2/hx1+2/hy1)^-1*( ( u(i-1,N+1) + u(i+1,N+1) )/hx1 + ( 2*u(i,N) )/hy1 - f(x(i),y(N)) );
    end
    error = max(max(abs(u-u2)));
    u2 = u;
    k = k+1;
end
figure(1);
surf(y,x,u2);
xlabel('x');
ylabel('y');
title('The approximate solution')
