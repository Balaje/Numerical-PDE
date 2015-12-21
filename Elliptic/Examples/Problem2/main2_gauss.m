% Program to solve Problem 4 of Elliptic PDE LabSheet

a = 0; 
b = 1;
c = 0;
d = 2;
N = 40;
M = 40;

hx = (b-a)/N;
hx1 = hx*hx;
hy = (d-c)/M;
hy1 = hy*hy;

x = a:hx:b;
y = c:hy:d;

mul = 1/(2/hx1+2/hy1);
u2 = zeros(N+1,M+1);
for i=1:N+1
    u2(i,1) = (x(i))^2;
    u2(i,M+1) = (x(i)-2)^2;
end
for i=1:M+1
    u2(1,i) = (y(i))^2;
    u2(N+1,i) = (y(i)-1)^2;
end

error = 1000;
tol = 10^-7;
k  = 0;

u = u2;
while error > tol
    for i=2:N
        for j=2:M
            u(i,j) = mul*( (u(i,j-1) + u(i,j+1))/hy1 +...
                (u(i-1,j)+u(i+1,j))/hx1 - 4); % f(x,y) = 4
        end
    end
    error = max(max(abs(u-u2)));
    u2 = u;
    k = k+1;
end

for i=1:M+1
    for j=1:N+1
        exact(i,j) = (x(i)-y(j))^2;
    end
end
figure(1);
surf(y,x,u2);
xlabel('x');
ylabel('y');
title('The approximate solution')

figure(2);
surf(y,x,exact);
xlabel('x');
ylabel('y');
title('The exact solution')
