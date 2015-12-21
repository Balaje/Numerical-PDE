function [x,exact] = genexact(eps,t0,tf,x0,xf)

N = 1000;
M = 1000;

h = (xf-x0)/N;
delt = (tf-t0)/M;

mu = delt/(2*h^2);
lamb = delt/(4*h);

x = zeros(N-1,1);
u = zeros(N-1,1);

for j=1:N-1
    x(j) = x0 + j*h;
    if((x(j)) < 1)
        u(j) = sin(pi*x(j));
    end
end

F = zeros(N-1,1);
U = ones(N-1,1);
J = zeros(N-1);
error = 100;
tol = 10^-5;

t = t0;
while error > tol
    for k=1:M
        F(1) = (1+2*mu*eps)*U(1) - mu*eps*U(2) + lamb*f(U(2)) - ...
            ( (1-2*mu*eps)*u(1) + eps*mu*u(2) - lamb*f(u(2)) );
        for j=2:N-2
            F(j) = -eps*mu*U(j-1) - lamb*f(U(j-1)) + (1+2*mu*eps)*U(j) - mu*eps*U(j+1) + lamb*f(U(j+1)) - ...
                ( lamb*f(u(j-1)) + eps*mu*u(j-1) + (1-2*mu*eps)*u(j) + eps*mu*u(j+1) - lamb*f(u(j+1)) );
        end
        F(N-1) = (1+2*mu*eps)*U(N-1) - mu*eps*U(N-2) + lamb*f(U(N-2)) - ...
            ( (1-2*mu*eps)*u(N-1) + eps*mu*u(N-2) - lamb*f(u(N-2)) );
        % Jacobian Matrix
        J(1,1) = 1+2*mu*eps;
        J(N-1,N-1) = 1+2*mu*eps;
        for j=2:N-2
            J(j,j) = 1+2*mu*eps;
            J(j-1,j) = -mu*eps + lamb*df(U(j));
            J(j,j-1) = -mu*eps - lamb*df(U(j-1));
        end
        J(N-1,N-2) = -mu*eps - lamb*df(U(N-2));
        J(N-2,N-1) = -mu*eps + lamb*df(U(N-1));
    
        DELT = J\F;
        U1 = U - DELT;
        error = max(abs(U1-U));
        u = U1;
        U = U1;
        t = t+delt;
        
    end
end
exact = U1;
end

