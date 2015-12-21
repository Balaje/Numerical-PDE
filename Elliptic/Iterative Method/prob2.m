% Jacobi Method for solving elliptic problem

a = 0; b = 1;
c = 0; d = 1;
N = 10;
hx = (b-a)/N;
M = 10;
hy = (d-c)/M;

x = a:hx:b;
y = c:hy:d;

u2 = zeros(N+1);
for j=1:N+1
    u2(j,1) = 0;
    u2(1,j) = 0;
    u2(j,N+1) = y(j);
    u2(N+1,j) = y(j);
end

error = 100;
tol = 10^-8;

u = u2;
k=0;
omega = 1.5;
while error > tol
    for i=2:N
        for j=2:N
            %u2(i,j) = 1/(2/hx^2+2/hy^2)*((u(i-1,j)+u(i+1,j))/hx^2 + (u(i,j-1)+u(i,j+1))/hy^2);
            u(i,j) = 1/(2/hx^2+2/hy^2)*((u(i-1,j)+u(i+1,j))/hx^2 + (u(i,j-1)+u(i,j+1))/hy^2);
            %u(i,j) = (1-omega)*u2(i,j) + omega*(1/(2/hx^2+2/hy^2)*((u(i-1,j)+u(i+1,j))/hx^2 + (u(i,j-1)+u(i,j+1))/hy^2));
        end
    end
    error = max(max(abs(u2-u)));
    k=k+1;
    %u = u2;
    u2 = u;
end

mesh(y,x,u2);