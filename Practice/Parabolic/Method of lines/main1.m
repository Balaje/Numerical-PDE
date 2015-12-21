% Neumann problem method of lines

global h N x

N = 100;

a = 0; b = 1;
h = (b-a)/N;
t0 = 0; tf = 1;
x = zeros(N,1);
u0 = zeros(N,1);

for j=1:N
    x(j) = a+j*h;
    u0(j) = sin(pi^2*x(j));
end

tspan = t0:0.01:tf;

[T U] = ode23s('mol1',tspan,u0);

exact = zeros(N,1);
for j=1:N
    exact(j) = exp(-tf)*sin(pi^2*x(j));
end

plot(x,U(size(T,1),:),x,exact);