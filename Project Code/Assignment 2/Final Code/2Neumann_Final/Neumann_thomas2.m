% Program to solve a Neumann elliptic problem using slow fourier transform
% and solving system of equations using the Thomas alogorithm
%Lab-Sheet Question 3. Change boundary conditions.

clear all
close all
clc

a = 0; b = 2;
c = 0; d = 1;

M = 300;

h = (b-a)/M;
N = fix((d-c)/h);

x = a+h:h:b;
y = c+h:h:d;

U = zeros(N,M);
A = zeros(N);
g = zeros(M,N);

tic
g(1,1) = f(x(1),y(1)) - 1/h^2 - 1/h^2;
for j=2:N-1
    g(1,j) = f(x(1),y(j)) - 1/h^2;
end
g(1,N) = f(x(1),y(N)) - 1/h^2 - (2*x(1)*exp(x(1)))/h;

g(M,1) = f(x(M),y(1)) - 1/h^2 - (2*y(1)*exp(2*y(1)))/h;
for j=2:N-1
    g(M,j) = f(x(M),y(j)) - (2*y(j)*exp(2*y(j)))/h;
end
g(M,N) = f(x(M),y(N)) - (2*x(M)*exp(x(M)))/h - (2*y(N)*exp(2*y(N)))/h;

for k=2:M-1
    g(k,1) = f(x(k),y(1)) - (1/h^2);
        for j=2:N-1
        g(k,j) = (f(x(k),y(j)));
        end
    g(k,N) = f(x(k),y(N)) - (2*x(k)*exp(x(k)))/h;
end

G = h^2*mydst(g);

% k=1:N-1
for k=1:M    
    A(1,1) = -(2+4*sin((2*k-1)*pi/(4*M))^2);
    A(1,2) = 1;
    A(N,N-1) = 2;
    A(N,N) = -(2+4*sin((2*k-1)*pi/(4*M))^2);
    for j=2:N-1
        A(j,j) = -(2+4*sin((2*k-1)*pi/(4*M))^2);
        A(j,j-1) = 1;
        A(j,j+1) = 1;
    end
    U(:,k) = thomasSolver(A,G(k,:)');
    %U(:,k) = A\(G(k,:)');
end

u = myidst(U');
toc

figure(1)
mesh(y,x,u);
title('Approximate solution using Fourier Method (Thomas Algorithm)');
xlabel('x');
ylabel('y');

for i=1:M
    for j=1:N
        usol(i,j)=exp(x(i)*y(j));
    end
end

figure(2);
mesh(y,x,usol);
title('The exact solution')

error=max(max(usol-u))