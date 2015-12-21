% Program to solve a dirichlet elliptic problem using iterative method 

clear
N = 100;
tol = 1e-6;
a = 0;
b = 1;
c = 0;
d = 1;
h = (b-a)/N;
h1 = h*h;
x = a:h:b;
y = c:h:d;
u2 = zeros(N+1,N+1);
u1 = zeros(N+1,N+1);
usol=zeros(N+1,N+1);

tic
for i = 1:N+1
    u2(1,i) = 0.1*cos(2*pi*y(i)); 
    u2(N+1,i) = 0.1*cos(2*pi*y(i)); 
    u2(i,1) = 0.1; 
    u2(i,N+1) = 0.1; 
end

k = 0;
while max(max(abs(u1-u2))) > tol
    k = k+1;
    u1 = u2;
    for i = 2:N
        for j = 2:N
            u2(i,j) = (u1(i-1,j) + u1(i+1,j) + u1(i,j-1) + u1(i,j+1)-h1*f(x(i),y(j)))/4;
        end
    end
end
toc

figure(3);
mesh(x,y,u2);
title('The approximate solution');
xlabel('x');
ylabel('y');

