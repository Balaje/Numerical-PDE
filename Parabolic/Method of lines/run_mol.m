% Program to implement the 
% Method of Lines
% Method to solve heat conduction equation

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %   u_t = alpha*u_xx + f(t,x) a<x<b, t>0 %
        %   u(x,0) = g(x)                        % 
        %   u(a,t) = 0,   u(b,t) = 0             %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       

clear all
close all

global h N x
N = 20;

a = 0; b = 1;
t0 = 0; tfinal = 1;

h = (b-a)/N;

x = zeros(N-1,1);
u0 = zeros(N-1,1);

for j=1:N-1
    x(j) = j*h;
    u0(j) = sin(pi*x(j));
end

tspan = t0:0.1:tfinal;

[t1 y] = ode15s('main_mol',tspan,u0);

plot(x,y(11,:));