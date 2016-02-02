% KDV Crank Nicolson with periodic boundary condition
% u(-20,t) = u(20,t)

N = 200;
M = 1000;

A = 10;
m = sqrt(A/2);

a = -20; b = 20;
t0 = 0; tf = 1;

h = (b-a)/N;
k = (tf-t0)/M;

mu = 6*k/(4*h);
tau = k/(4*h^3);

x = a:h:b;

u0 = A*(sech(m*x)).^2;

F = zeros(N+1,1);
J = zeros(N+1);
error = 100;
tol = 10^-6;

U = ones(N+1,1);
while error > tol
    for p=1:M
        F(1) = U(1) - u0(1) + mu*(f(U(2))-f(U(N)) + f(u0(2))-f(u0(N))) + ...
            tau*(U(3)-2*U(2)+2*U(N)-U(N-1) + u0(3)-2*u0(2)+2*u0(N)-u0(N-1));    
        
        F(2) = U(2) - u0(2) + mu*(f(U(3))-f(U(1)) + f(u0(3))-f(u0(1))) + ...
            tau*(U(4)-2*U(3)+2*U(1)-U(N) + u0(4)-2*u0(3)+2*u0(1)-u0(N));
        
        F(N) = U(N) - u0(N) + mu*(f(U(N+1))-f(U(N-1)) + f(u0(N+1))-f(u0(N-1))) + ...
            tau*(U(2)-2*U(N+1)+2*U(N-1)-U(N-2) + u0(2)-2*u0(N+1)+2*u0(N-1)-u0(N-2));
        
        F(N+1) = U(N+1) - u0(N+1) + mu*(f(U(2))-f(U(N)) + f(u0(2))-f(u0(N))) + ...
            tau*(U(3)-2*U(2)+2*U(N)-U(N-1) + u0(3)-2*u0(2)+2*u0(N)-u0(N-1));
        
        for j=3:N-1
            F(j) = U(j) - u0(j) + mu*(f(U(j+1))-f(U(j-1)) + f(u0(j+1))-f(u0(j-1))) + ...
                tau*(U(j+2)-2*U(j+1)+2*U(j-1)-U(j-2) + u0(j+2)-2*u0(j+1)+2*u0(j-1)-u0(j-2));
        end
        
        % Jacobian
        % Jacobian
        J(1,1) = 1;
        J(1,2) = mu*df(U(2)) - 2*tau;
        J(1,3) = tau;
        J(1,N-1) = -tau;
        J(1,N) = -mu*df(U(N))+2*tau;
        
        J(2,1) = -mu*df(U(1)) + 2*tau;
        J(2,2) = 1;
        J(2,3) = mu*df(U(3)) - 2*tau;
        J(2,4) = tau;
        J(2,N) = -tau;
        
        for j=3:N-1
            J(j,j) = 1;
            J(j,j+2) = tau;
            J(j,j-2) = -tau;
            J(j-1,j) = mu*df(U(j))-2*tau;
            J(j,j-1) = -mu*df(U(j-1))+2*tau;
        end
        J(N,2) = tau;
        J(N,N-1) = -mu*df(U(N-1))+2*tau;
        J(N,N) = 1;
        J(N,N-2) = -tau;
        J(N,N+1) = mu*df(U(N)) - 2*tau;
        
        J(N+1,2) = mu*df(U(2))-2*tau;
        J(N+1,3) = tau;
        J(N+1,N-1) = -tau;
        J(N+1,N) = -mu*df(U(N))+2*tau;
        J(N+1,N+1) = 1;
      
        DELT = J\F;
        U2 = U - DELT;
        error = max(abs(U2-U));
        U = U2;
        u0 = U2;
        plot(x,U)
        axis([a,b,-2,10])
        pause(0.001)
    end
end