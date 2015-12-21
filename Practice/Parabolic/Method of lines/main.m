% Program to solve a parabolic problem using method of lines

global h N x

N = 200;

a = 0; b = 2;
t0 = 0; tf = 1;

h = (b-a)/N;

x = zeros(N-1,1);
u0 = zeros(N-1,1);

for j=1:N-1
    x(j) = j*h;
    u0(j) = sin((pi/2)*x(j));
end

tspan = t0:0.01:tf;

[T U] = ode15s('mol',tspan,u0);

exact = zeros(N-1,1);
for j=1:N-1
    exact(j) = exp(-pi^2/4*tf)*sin(pi/2*x(j));
end

plot(x,U(size(T,1),:),x,exact)