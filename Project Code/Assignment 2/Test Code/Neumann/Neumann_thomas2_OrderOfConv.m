% Program to solve a Neumann elliptic problem using slow fourier transform
% with order of convergenece 
% and ny solving system of equations using the Thomas alogorithm
%Lab-Sheet Question 3. Change boundary conditions.

clear

a = 0; b = 2;
c = 0; d = 1;

M = 10;
tic
for p=1:5
    h(p) = (b-a)/M;
    N = fix((d-c)/h(p));

    x = a+h(p):h(p):b;
    y = c+h(p):h(p):d;

    U = zeros(N,M);
    A = zeros(N);
    g = zeros(M,N);

   
    g(1,1) = f(x(1),y(1)) - 1/h(p)^2 - 1/h(p)^2;
    for j=2:N-1
        g(1,j) = f(x(1),y(j)) - 1/h(p)^2;
    end
    g(1,N) = f(x(1),y(N)) - 1/h(p)^2 - (2*x(1)*exp(x(1)))/h(p);

    g(M,1) = f(x(M),y(1)) - 1/h(p)^2 - (2*y(1)*exp(2*y(1)))/h(p);
    for j=2:N-1
        g(M,j) = f(x(M),y(j)) - (2*y(j)*exp(2*y(j)))/h(p);
    end
    g(M,N) = f(x(M),y(N)) - (2*x(M)*exp(x(M)))/h(p) - (2*y(N)*exp(2*y(N)))/h(p);

    for k=2:M-1
        g(k,1) = f(x(k),y(1)) - (1/h(p)^2);
            for j=2:N-1
            g(k,j) = (f(x(k),y(j)));
            end
        g(k,N) = f(x(k),y(N)) - (2*x(k)*exp(x(k)))/h(p);
    end

    G = h(p)^2*mydst(g);

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

    for i=1:M
        for j=1:N
            usol(i,j)=exp(x(i)*y(j)); %exact solution
        end
    end

    error(p) = max(max(abs(usol-u)));

    M=2*M;

end
error

for j=1:p-1
    order(j)= log(error(j)/error(j+1))/log(2);
end
order

toc

figure(1)
surf(y,x,u);
title('Approximate solution using Fourier Method (Thomas Algorithm)');
xlabel('x');
ylabel('y');


figure(2);
surf(y,x,usol);
title('The exact solution')
