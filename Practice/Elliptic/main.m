% Program to solve an elliptic problem using Jacobi Iteration 

N = 20;
a = 0; b = 1;
c = 0; d = 1;

h = (b-a)/N;

x = a:h:b;
y = c:h:d;

u0 = zeros(N+1);
for i=1:N+1
    u0(i,1) = sin(pi*x(i));
    u0(i,N+1) = -sin(pi*x(i));
    u0(1,i) = 0;
    u0(N+1,i) = 0;
end

tol = 10^-8;
error = 100;
u = u0;

while error > tol
    u0 = u;
    for i=2:N
        for j=2:N
            u(i,j) = (1/ ( 2/h^2 + 2/h^2*p(x(i),y(j)) - r(x(i),y(j)) ))*...
                ( (u0(i-1,j) + u0(i+1,j) )/h^2 + p(x(i),y(j))*(u0(i,j-1)+u0(i,j+1))/h^2 -...
                f(x(i),y(j)) );
        end
    end
    error = max(max(abs(u-u0)));
end

ue = zeros(N+1);
for i=1:N+1
    for j=1:N+1
        ue(i,j) = sin(pi*x(i))*cos(pi*y(j));
    end
end

figure(1);
mesh(y,x,u);
figure(2);
mesh(y,x,ue);