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

clc

% Input variables
alpha = 1/16;
N = 15; % Number of subintervals
k = 0.02; % Time step size
a = 0; b = 1; % Domain
t0 = 0; tf = 1;

% Calculate step size
h = (b-a)/N;

% Unconditionally stable
mu = alpha*k/h^2;

% Constructing the system of equations
% Declare a sparse matrix
A = sparse(N-1,N-1);
for i=1:N-1
    A(i,i) = 1+2*mu;  % Main diagonal
end
for i=1:N-2
    A(i+1,i) = -mu; % Lower diagonal
    A(i,i+1) = -mu; % Upper diagonal
end
% Right Hand Side
B = zeros(N-1,1);

M = fix((tf-t0)/k); %Total number of time steps
tl = t0:k:tf;
u0 = zeros(1,N-1);
x = zeros(1,N-1);

% Initial Conditions
for i=1:N-1
    x(i) = i*h;
    u0(i) = 2*sin(2*pi*x(i));
end

% Construct the system of equations
for i=1:M-1
    B(1) = u0(1) + k*f(tl(i+1),x(1));
    B(N-1) = u0(N-1)+ k*f(tl(i+1),x(N-1));
    
    for j=2:N-2
        B(j) = u0(j) + k*f(tl(i+1),x(j));
    end
    
    U_new = A\B;
    u0 = U_new;
end

% Compare with the exact solution
exact = zeros(1,size(u0,1));
for j=1:N-1
    exact(j) = 2*exp(-pi^2/4*tl(i))*sin(2*pi*x(j));
end
plot(x,u0,'b-',x,exact,'r-+');