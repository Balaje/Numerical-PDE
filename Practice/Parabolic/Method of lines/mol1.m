function dU = mol1(t,u)

global h N x

dU = zeros(N,1);

dU(1) = (-2*u(1)+u(2))/h^2 + f(x(1),t);
for j=2:N-1
    dU(j) = (u(j-1) - 2*u(j) + u(j+1))/h^2 + f(x(j),t);
end
dU(N) = (2*u(N-1)+ 2*h*exp(-t)*pi^2*cos(pi^2) - 2*u(N))/h^2 + f(x(N),t);