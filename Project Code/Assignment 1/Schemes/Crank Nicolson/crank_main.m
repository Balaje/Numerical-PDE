% Program to solve Burger's Equation using Crank Nicolson Method

clear all
eps = 0.09;

x0 = 0; xf = 1;
t0 = 0; tf = 1;

[x1,exact] = genexact(eps,t0,tf,x0,xf);

N = 10;
M = 1000;
z = 4;

error2 = zeros(z,1);
order = zeros(z-1,1);
h = zeros(z,1);

for p=1:z
    h(p) = (xf-x0)/N;
    delt = (tf-t0)/M;

    mu = delt/(2*h(p)^2);
    lamb = delt/(4*h(p));

    x = zeros(N-1,1);
    u = zeros(N-1,1);

    for j=1:N-1
        x(j) = x0 + j*h(p);
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
            for j=1:N-1
                J(j,j) = 1+2*mu*eps;
            end
            for j=2:N-1
                J(j-1,j) = -mu*eps + lamb*df(U(j));
                J(j,j-1) = -mu*eps - lamb*df(U(j-1));
            end
            DELT = J\F;
            U1 = U - DELT;
            error = max(abs(U1-U));
            u = U1;
            U = U1;
            t = t+delt;
            %plot(x,U1)
             %str = ['At time t=',num2str(t)];
             %title(str);
             %axis([x0,xf,0,1])
             %pause(0.01);
        end
    end
    
    error1 = zeros(gcd(size(x1,1)+1,size(x,1)+1),1);
    for i=1:N-1
        for j=1:size(exact,1)-1
            if(x(i) == x1(j))
                error1(i) = abs(exact(j)-U1(i));
            end
        end
    end
    error2(p) = max(error1);
    N = N*2;
end

for j=1:p-1
    order(j) = log(error2(j)/error2(j+1))/log(2);
end