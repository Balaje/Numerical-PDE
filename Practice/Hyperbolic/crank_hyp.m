% Program to implement the Crank Nicolson Scheme for a hyperbolic equation.
%           u_t + s*u_x = 0;          0<x<4, t>0
%           u(x,0) = u0(x) = 0  , 0<x<2
%                          = sin(pi*(x-2)), 2<x<4

% Wave speed
s = 1;

a = 0; b = 2;
N = 10;
h = (b-a)/N;

t0 = 0; tf = 1;
M = 100;
k = (tf-t0)/M;

x = a+h:h:b-h;
u0 = zeros(N-1,1);

for j=1:N-1
   if(x(j) < 0)
       
   end
end

mu = -(s*k)/(4*h);

A = zeros(N-1);
B = A;
A(1,1) = 1;
B(1,1) = 1;
for j=2:N-1
    A(j,j) = 1;
    A(j,j-1) = mu;
    A(j-1,j) = -mu;
    B(j,j) = 1;
    B(j,j-1) = -mu;
    B(j-1,j) = mu;
end

u = zeros(N-1,1);
for j=1:M
    u = A\(B*u0);
    u0 = u;
end

exact = zeros(N-1,1);
for j=1:N-1
   if(x(j)-tf < 2)
       exact(j) = 0;
   elseif(2 <= x(j)-tf && x(j)-tf < 4)
       exact(j) = sin(pi*(x(j)-tf));
   end
end

plot(x,exact,x,u);