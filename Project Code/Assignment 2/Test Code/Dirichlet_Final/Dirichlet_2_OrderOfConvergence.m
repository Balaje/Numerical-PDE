% Program to solve a dirichlet elliptic problem using fft and 
% solving system of equations using the Thomas alogorithm

clear all
close all
clc
format long

a = 0; b = 2;
c = 0; d = 1;

M = 10;

% Start timer
tic
for p=1:5
    h(p) = (b-a)/M;

    x = a+h(p):h(p):b-h(p); %Assuming equal grid spacing in x and y direction
    y = c+h(p):h(p):d-h(p);

    N=fix((d-c)/h(p));

    U = zeros(N-1,M-1);
    A = zeros(M-1,N-1);
    g = zeros(M-1,N-1);


    % Computing the rhs vectors
    g(1,1) = (f(x(1),y(1)) - 1/h(p)^2 - 1/h(p)^2);
    for j=2:N-2
        g(1,j) = (f(x(1),y(j)) - 1/h(p)^2);
    end
    g(1,N-1) = (f(x(1),y(N-1)) - 1/h(p)^2 - exp(x(1))/h(p)^2);

    g(M-1,1) = (f(x(M-1),y(1)) - exp(2*y(1))/h(p)^2 - 1/h(p)^2);
    for j=2:N-2
        g(M-1,j) = (f(x(M-1),y(j)) - exp(2*y(j))/h(p)^2);
    end
    g(M-1,N-1) = (f(x(M-1),y(N-1)) - exp(2*y(N-1))/h(p)^2 - exp(x(M-1))/h(p)^2);

    for i=2:M-2
        g(i,1) = (f(x(i),y(1)) - 1/h(p)^2);
            for j=2:N-2
                g(i,j) = (f(x(i),y(j)));
            end
        g(i,N-1) = (f(x(i),y(N-1)) - exp(x(i))/h(p)^2);
    end

    % Computing the discrete sine transform
    %----------------
    G=h(p)^2*mydst(g);
    %G = h(p)^2*dst(g);
    %----------------

    % k = 1;
    A(1,1) = -(2+4*sin(pi/(2*M))^2);
    for j=2:N-1
        A(j,j) = -(2+4*sin(pi/(2*M))^2);
        A(j,j-1) = 1;
        A(j-1,j) = 1;
    end
    U(:,1) = thomasSolver(A,G(1,:)');

    % k=2:M-2
    for k=2:M-2    
        A(1,1) = -(2+4*sin(k*pi/(2*M))^2);
        for j=2:N-1
            A(j,j) = -(2+4*sin(k*pi/(2*M))^2);
            A(j,j-1) = 1;
            A(j-1,j) = 1;
        end
        U(:,k) = thomasSolver(A,G(k,:)');
    end

    % Computing the tridiagonal [A] matrix and the fourier coefficients
    %k=M-1
    A(1,1) = -(2+4*sin((M-1)*pi/(2*M))^2);
    for j=2:N-1
        A(j,j) = -(2+4*sin((M-1)*pi/(2*M))^2);
        A(j,j-1) = 1;
        A(j-1,j) = 1;
    end
    U(:,M-1) = thomasSolver(A,G(M-1,:)');


    u=zeros(M-1,N-1);
    % Computing the inverse discrete sine transform
    %------------------
    u = myidst(U');
    %u=idst(U');
    %------------------

    usol=zeros(M-1,N-1);
    for i=1:M-1
        for j=1:N-1
            usol(i,j)=exp(x(i)*y(j));
        end
    end



    M=2*M;
    
    error(p) = max(max(abs(u-usol)));
    
end
error

for j=1:p-1
    order(j) = log(error(j)/error(j+1))/log(2);
end

order

toc

    figure(1)
    mesh(y,x,u);
    title('Approximate solution using Fourier Method (Thomas Algorithm)');
    xlabel('x');
    ylabel('u');
    
    figure(2)
    mesh(y,x,usol);
    title('The exact solution');
    xlabel('x');
    ylabel('u');
    
    figure(3)
    plot(log(h),log(error));
    title('log(error) vs log(h)');
    xlabel('log(h)');
    ylabel('log(error)');