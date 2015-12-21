% Program to implement the 
% Backward in time-Central in space
% Method to solve heat conduction equation
        
clear all 
close all

format long

clc

% Input variables
alpha = 1;
N = 10; % Number of subintervals

error = zeros(10,1);
order = zeros(9,1);
h = zeros(10,1);
% Calculate step size
for p=1:1
    M = 1000; % Time step size
    a = 0; b = 1; % Domain
    t0 = 0; tf = 1;
    h(p) = (b-a)/N;
    k = (tf-t0)/M;
    % Unconditionally stable
    mu = alpha*k/(2*h(p)^2);

    % Constructing the system of equations
    A = zeros(N+1,N+1);
    C = zeros(N+1,N+1);

    A(1,1) = -2*mu*h(p)+1+2*mu;
    A(1,2) = -2*mu;
    A(N+1,N) = -2*mu;

    C(1,1) = 1+2*mu*h(p)-2*mu;
    C(1,2) = 2*mu;
    C(N+1,N) = 2*mu;

    for i=2:N+1
        A(i,i) = 1+2*mu;  % Main diagonal
        C(i,i) = 1-2*mu;
    end
    for i=1:N-1
        A(i+1,i) = -mu; % Lower diagonal
        A(i+1,i+2) = -mu; % Upper diagonal
        C(i+1,i) = mu;
        C(i+1,i+2) = mu;
    end

    % Right Hand Side
    b = zeros(N+1,1);

    u0 = zeros(N+1,1);
    x = zeros(N+1,1);

    % Initial Conditions
    for i=1:N+1
        x(i) = (i-1)*h(p);
        u0(i) = x(i)*(1-x(i));
    end

    % Construct the system of equations
    t=t0;
    for i=1:M-1
        b(1) = 0.5*k*( f(t+k,x(1)) + f(t,x(1)) ) - 2*mu*h(p)*exp(-(t+k)) - 2*mu*h(p)*exp(-t);
        b(N+1) = 0.5*k*(f(t+k,x(N+1))+f(t,x(N+1))) - 2*h(p)*mu*exp(-(t+k)) - 2*mu*h(p)*exp(-t);
    
        for j=2:N
            b(j) = 0.5*k*(f(t+k,x(j))+f(t,x(j)));
        end
    
        U_new = A\(C*u0+b);
        u0 = U_new;
        t=t+k;
    end

    % Compare with the exact solution
    exact = zeros(size(u0,1),1);
    for j=1:N+1
        exact(j) = exp(-t)*x(j)*(1-x(j));
    end

    error(p) = norm((u0-exact),Inf);
    N = N*2;
    
    for j=1:p-1
        order(j) = log(error(j)/error(j+1))/log(2);
    end

end

figure(1);
plot(x,u0,'b-+',x,exact,'r-');
figure(2);
plot(log(h),log(error));