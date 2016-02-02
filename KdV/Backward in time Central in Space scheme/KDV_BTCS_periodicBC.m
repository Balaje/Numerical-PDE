% KDV Backward in time with periodic boundary condition
% u(-20,t) = u(20,t)

N = 400;
M = 10000;

A = 8;
m = sqrt(A/2);

a = -20; b = 20;
t0 = 0; tf = 1;

h = (b-a)/N;
k = (tf-t0)/M;

mu = 3*k/h;
tau = k/h^3;

x = a:h:b;

U0 = A*(sech(m*x)).^2;

F = zeros(N+1,1);
J = zeros(N+1,N+1);
error = 100;
tol = 10^-6;

U=ones(N+1,1);
while error>tol
    for j=1:M

        F(1) = U(1) - U0(1) + mu*[f(U(2)) - f(U(N))] + tau*[-0.5*U(N-1) + U(N)...
            - U(2) + 0.5*U(3)];

        F(2) = U(2) - U0(2) + mu*[f(U(3)) - f(U(1))] + tau*[-0.5*U(N) + U(1)...
            - U(3) + 0.5*U(4)];
        
        for i=3:N-1
            F(i) = U(i) - U0(i) + mu*[f(U(i+1)) - f(U(i-1))] + tau*[-0.5*U(i-2) + U(i-1)...
                - U(i+1) + 0.5*U(i+2)];
        end

        F(N) = U(N) + U0(N) + mu*[f(U(N+1)) - f(U(N-1))] + tau*[-0.5*U(N-2)...
             + U(N-1) - U(N) + 0.5*U(2)];

        F(N+1) = U(N+1) - U0(N+1) + mu*[f(U(2)) - f(U(N))] + tau*[-0.5*U(N-1)...
            + U(N) - U(2) + 0.5*U(3)];


        % Jacobian
        % Jacobian
        J(1,1) = 1;
        J(1,2) = mu*fp(U(2)) - tau;
        J(1,3) = tau/2;
        J(1,N-1) = -tau/2;
        J(1,N) = tau;
        
        J(2,1) = -mu*fp(U(1)) + tau;
        J(2,2) = 1;
        J(2,3) = mu*fp(U(3)) - tau;
        J(2,4) = tau/2;
        J(2,N) = -tau/2;
        
        for j=3:N-1
            J(j-2,j) = -tau/2;
            J(j,j-1) = -mu*fp(U(j-1)) + tau;
            J(j,j) = 1;
            J(j,j+1) = mu*fp(U(j+1)) - tau;
            J(j,j+2) = tau/2;
        end
        
        J(N,2) = tau/2;
        J(N,N-2) = -tau/2; 
        J(N,N-1) = -mu*fp(U(N-1)) + tau;
        J(N,N) = 1;
        J(N,N+1) = mu*fp(U(N+1)) - tau;
        
        J(N+1,2) = mu*fp(U(2)) - tau;
        J(N+1,3) = tau/2;
        J(N+1,N-1) = -tau/2;
        J(N+1,N) = -mu*fp(U(N)) + tau;
        J(N+1,N+1) = 1;
      
        DELT = J\F;
        U1 = U - DELT;
        error = max(abs(U1-U));
        U = U1;
        U0 = U1;
        
        plot(x,U);
        axis([a,b,-2,10])
        pause(1E-10);
        
    end
end
