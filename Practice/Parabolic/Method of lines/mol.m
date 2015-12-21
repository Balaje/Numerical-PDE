function dU = mol(t,u)

global h N x

dU = zeros(N-1,1);

dU(1) = (-2*u(1) + u(2))/h^2;

for j=2:N-2
    dU(j) = (u(j-1) - 2*u(j) + u(j+1))/h^2;
end

dU(N-1) = (u(N-2) - 2*u(N-1))/h^2;

