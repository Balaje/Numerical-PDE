% Program to implement the 
% Backward in time-Central in space
% Method to solve heat conduction equation

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %   u_t = alpha*u_xx + f(t,x) a<x<b, t>0 %
        %   u(x,0) = g(x)                        % 
        %   u(a,t) = 0,   u(b,t) = 0             %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
clear all 
close all

% Input variables
alpha = 1;
N = 10; % Number of subintervals
M = 10; % Time step size
a = 0; b = 1; % Domain
t0 = 0; tf = 1;

% Calculate step size
h = (b-a)/N;
k = (tf-t0)/M;
% Unconditionally stable
mu = alpha*k/(2*h^2);

% Constructing the system of equations
% Declare a sparse matrix
A = sparse(N-1,N-1);
C = sparse(N-1,N-1);
for i=1:N-1
    A(i,i) = 1+2*mu;  % Main diagonal
    C(i,i) = 1-2*mu;
end
for i=1:N-2
    A(i+1,i) = -mu; % Lower diagonal
    A(i,i+1) = -mu; % Upper diagonal
    C(i+1,i) = mu;
    C(i,i+1) = mu;
end

% Right Hand Side
b = zeros(N-1,1);

u0 = zeros(N-1,1);
x = zeros(1,N-1);

% Initial Conditions
for i=1:N-1
    x(i) = i*h;
    u0(i) = sin(pi*x(i));
end

% Construct the system of equations
t=t0;
for i=1:M-1
    %b(1) = k*(f(t+0.5*k,x(1)));
    %b(N-1) = k*(f(t+0.5*k,x(N-1)));
    
    for j=1:N-1
        b(j) = k*(f(t+0.5*k,x(j)));
    end
    
    U_new = A\(C*u0+b);
    u0 = U_new;
    t=t+k;
end

% Compare with the exact solution
exact = zeros(1,size(u0,1));
for j=1:N-1
    exact(j) = exp(-t)*sin(pi*x(j));
end
plot(x,u0,'b-',x,exact,'r-+');