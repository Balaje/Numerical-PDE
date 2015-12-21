% Problem u_xx+u_yy = -2*pi^2*u
% u = 0 on boundary.

N = 50;
a = 0; b = 1;
c = 0; d = 1;
h = (b-a)/N;

x = a:h:b;
y = c:h:d;

u0 = zeros(N+1);
for j=1:N+1
    u0(j,1) = sin(pi*x(j));
    u0(j,N+1) = -sin(pi*x(j));
    u0(1,j) = 0;
    u0(N+1,j) = 0;
end
u = u0;

error = 100;
tol = 10^-8;

while error > tol
    for i=2:N
        for j=2:N
            u(i,j) =h^2/4*( (u0(i-1,j)+u0(i+1,j)+u0(i,j-1)+u0(i,j+1))/h^2 - f1(x(i),y(j)) );
        end
    end
    error = max(max(abs(u-u0)));
    u0 = u;
end

figure(1)
surf(y,x,u);
title('Approx');

ue = zeros(N+1);
for i=1:N+1
    for j=1:N+1
        ue(i,j) = sin(pi*x(i))*cos(pi*y(j));
    end
end
figure(2)
surf(y,x,ue);
title('Exact');