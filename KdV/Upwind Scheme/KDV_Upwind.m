% First order upwind scheme for the KdV equation
% Periodic boundary condition

A = 2;
m = sqrt(A/2);
a = -20; b = 20;
h = 0.08;
N = fix((b-a)/h);

tf = 1; t0 = 0;
k = 0.000001;
M = fix((tf-t0)/k);

x = a:h:b;

u0 = A*sech(m*x).^2;

U = zeros(N+1,1);

mu = 6*k/h;
tau = k/(2*h^3);
for p=1:M
    U(1) = u0(1)-mu*u0(1)*(f1(u0(1),u0(2))-f1(u0(N),u0(1))) + tau*(u0(3)-2*u0(2)+2*u0(N)-u0(N-1));
    U(2) = u0(2)-mu*u0(2)*(f1(u0(2),u0(3))-f1(u0(1),u0(2))) + tau*(u0(4)-2*u0(3)+2*u0(1)-u0(N));
    for j=3:N-1
        U(j) = u0(j) - mu*u0(j)*(f1(u0(j),u0(j+1))-f1(u0(j-1),u0(j))) + tau*(u0(j+2)-2*u0(j+1)+2*u0(j-1)-u0(j-2));
    end
    U(N) = u0(N)-mu*u0(N)*(f1(u0(N),u0(N+1))-f1(u0(N-1),u0(N))) + tau*(u0(2)-2*u0(N+1)+2*u0(N-1)-u0(N-2));
    U(N+1) = u0(N+1)-mu*u0(N+1)*(f1(u0(N+1),u0(2))-f1(u0(N),u0(N+1))) + tau*(u0(3)-2*u0(2)+2*u0(N)-u0(N-1));
    
    plot(x,U,'b');
    axis([a,b,-A,A])
    pause(0.001);
    u0 = U;
end