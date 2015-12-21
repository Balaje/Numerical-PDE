% Solve Burger's Equaiton using backward euler

clear all

a = 0; b = 1;
t0 = 0; tf = 1;
eps = 0.09;

N = 10;
M = 1000;

[x1,exact] = genexact(a,b,t0,tf,eps);

z = 6;
order = zeros(z-1,1);
error2 = zeros(z,1);
h = zeros(z,1);

for p=1:z
    h(p) = (b-a)/N;
    delt = (tf-t0)/M;

    lambda = delt/(2*h(p));
    mu = delt/(2*h(p)^2);

    x = zeros(N-1,1);
    u0 = zeros(N-1,1);

    for j=1:N-1
        x(j) = a+j*h(p);
        u0(j) = sin(pi*x(j));
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
%             if(p==2)
%                 plot(x,U1)
%                 axis([a,b,0,1.5]);
%                 pause(0.01);
%             end
            t = t+delt;
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
figure(1)
str = ['At time t=',num2str(t)];
plot(x,U1,x1,exact);
title(str)
xlabel('x');
ylabel('u(x,t)');
figure(2)
plot(log(h),log(error2));
