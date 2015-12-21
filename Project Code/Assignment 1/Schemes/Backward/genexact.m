function [x1,exact] = genexact(a,b,t0,tf,eps)
% Solve Burger's Equaiton using backward euler

N = 1000;
M = 1000;

h = (b-a)/N;
delt = (tf-t0)/M;

lambda = delt/(2*h);
mu = delt/(2*h^2);

x1 = zeros(N-1,1);
u0 = zeros(N-1,1);

for j=1:N-1
    x1(j) = a+j*h;
    u0(j) = sin(pi*x1(j));
end

t = t0;
U = ones(N-1,1);
F = zeros(N-1,1);
J = zeros(N-1);

error = 100;
tol = 10^-5;

while error > tol
    for k=1:M
        F(1) = U(1) + lambda*(f(U(2))) - mu*eps*(-2*U(1)+U(2)) - u0(1);
        for j=2:N-2
            F(j) = U(j) + lambda*(f(U(j+1)) - f(U(j-1))) - mu*eps*(U(j-1)-2*U(j)+U(j+1)) - u0(j);
        end
        F(N-1) = U(N-1) - lambda*(f(U(N-2))) - mu*eps*(-2*U(N-1)+U(N-2)) - u0(N-1);
        
        % Jacobian
        for j=1:N-1
            J(j,j) = 1+2*mu*eps;
        end
        for j=2:N-1
            J(j,j-1) = -lambda*df(U(j-1)) - mu*eps;
            J(j-1,j) = lambda*df(U(j)) - mu*eps;
        end
        
        DELTA = J\F;
        U1 = U - DELTA;
        error = max(abs(U1-U));
        U = U1;
        u0 = U1;
        t = t+delt;
    end
end

exact = U1;
end