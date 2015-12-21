% Program to solve a dirichlet problem using fft and 
% solving system of equations using the Gaussian Elimination
clear

a = 0; b = 1;
c = 0; d = 1;

% Assuming equal partition in x and y directions
N = 1000;

h = (b-a)/N;

x = a+h:h:b-h;
y = c+h:h:d-h;

U = zeros(N-1);
A = zeros(N-1);
g = zeros(N-1,N-1);

% Start timer
tic
% Computing the rhs vectors
g(1,1) = (f(x(1),y(1)) - 0.1*cos(2*pi*y(1))/h^2 - 0.1/h^2);
g(1,N-1) = (f(x(1),y(N-1)) - 0.1*cos(2*pi*y(N-1))/h^2 - 0.1/h^2);
for j=2:N-2
    g(1,j) = (f(x(1),y(j)) - 0.1*cos(2*pi*y(j))/h^2);
end

g(N-1,1) = (f(x(N-1),y(1)) - 0.1*cos(2*pi*y(1))/h^2 - 0.1/h^2);
g(N-1,N-1) = (f(x(N-1),y(N-1)) - 0.1*cos(2*pi*y(N-1))/h^2 - 0.1/h^2);
for j=2:N-2
    g(N-1,j) = (f(x(N-1),y(j)) - 0.1*cos(2*pi*y(j))/h^2);
end

for k=2:N-2
    g(k,1) = (f(x(k),y(1)) - 0.1/h^2);
    g(k,N-1) = (f(x(k),y(N-1)) - 0.1/h^2);
    for j=2:N-2
        g(k,j) = (f(x(k),y(j)));
    end 
end

% Computing the discrete sine transform
%----------------
G = h^2*dst(g);
%----------------


% Computing the [A] matrix and the fourier coefficients
% k = 1;
A(1,1) = -(2+4*sin(pi/(2*N))^2);
for j=2:N-1
    A(j,j) = -(2+4*sin(pi/(2*N))^2);
    A(j,j-1) = 1;
    A(j-1,j) = 1;
end
U(:,1) = A\G(1,:)';

% k=2:N-2
for k=2:N-2    
    A(1,1) = -(2+4*sin(k*pi/(2*N))^2);
    for j=2:N-1
        A(j,j) = -(2+4*sin(k*pi/(2*N))^2);
        A(j,j-1) = 1;
        A(j-1,j) = 1;
    end
    U(:,k) = A\G(k,:)';
end

%k=N-1
A(1,1) = -(2+4*sin((N-1)*pi/(2*N))^2);
for j=2:N-1
    A(j,j) = -(2+4*sin((N-1)*pi/(2*N))^2);
    A(j,j-1) = 1;
    A(j-1,j) = 1;
end
U(:,N-1) = A\G(N-1,:)';

% Computing the inverse discrete sine transform
%------------------
u = idst(U');
%------------------

% End timer
toc

figure(1)
mesh(x,y,u);
title('Approximate solution using Fourier Method');
xlabel('x');
ylabel('y');