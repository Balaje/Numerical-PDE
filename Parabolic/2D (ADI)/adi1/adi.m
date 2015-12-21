% Program to implement the ADI method for two dimensional
% Parabolic problem

% u_t = u_xx + u_yy | (x,y) in (0,1) x (0,1)

% u(x,y,0) = sin(pi*x)*sin(pi*y)
% u = 0 on the boundary

clear all
close all

format long
a1 = 0; % x-direction
b1 = 1;

c = 0; % y-direction
d = 1;

t0 = 0; tf = 0.0001; % Time
z = 3;

h = zeros(z,1);
N = 10; % Space Discretization [For non uniform discretization declare two
        %                       Space Variables]
M = 10; % Time discretization

% Loop to find order of convergence
for l = 1:z
    h(l) = (b1-a1)/N; % Discretization in the x direction. Discretization in the y direction (Assuming h_x = h_y)

    del_t = (tf-t0)/M;

    u0 = zeros(N-1,N-1);
    unew = u0;
    x = a1+h(l):h(l):b1-h(l); % Using the fact that discretization is same in both directions
    y = c+h(l):h(l):d-h(l);

    for j=1:N-1
        for k=1:N-1
            u0(j,k) = sin(pi*x(j))*sin(pi*y(k));
        end
    end

    % Defining the A matrix
    mu = del_t/(2*h(l)^2);
    A = sparse(N-1,N-1);
    for j=1:N-1
        A(j,j) = 1+2*mu;
    end
    for j=2:N-1
        A(j,j-1) = -mu;
        A(j-1,j) = -mu;
    end

    t = t0;

    for p=1:M
        b = zeros(N-1,1);
        %First part : Use 'unew' to store current data
        for i=1:N-1
            b(i) = (1-2*mu)*u0(i,1)+ mu*u0(i,2);
        end
        unew(1,:) = A\b;
    
        for j=2:N-2
            for i=1:N-1
                b(i) = mu*u0(i,j-1)+(1-2*mu)*u0(i,j)+mu*u0(i,j+1);
            end
            unew(j,:) = A\b;
        end
    
        for i=1:N-1
            b(i) = (1-2*mu)*u0(i,N-1)+mu*u0(i,N-2);
        end
        unew(N-1,:) = A\b; 
    
        %Second part : Update u0 and solve the second system
        for i=1:N-1
            b(i) = (1-2*mu)*unew(i,1)+mu*unew(i,2);
        end
        u0(1,:) = A\b;
    
        for j=2:N-2
            for i=1:N-1
                b(i) = mu*unew(i,j-1)+(1-2*mu)*unew(i,j)+mu*unew(i,j+1);
            end
            u0(j,:) = A\b;
        end
    
        for i=1:N-1
            b(i) = (1-2*mu)*unew(i,N-1)+mu*unew(i,N-2);
        end
        u0(N-1,:) = A\b;
    
        t = t+del_t;
    end

    exact = zeros(N-1,N-1);
    for i=1:N-1
        for j=1:N-1
            exact(i,j) = exp(-2*pi^2*t)*sin(pi*x(i))*sin(pi*x(j));
        end
    end

    error(l) = max(max(abs(u0-exact)));

    for j=1:l-1
        order(j) = log(error(j)/error(j+1))/log(2);
    end

    N=2*N;
end

figure(1)
surf(x,y,u0);
title('Approximate solution');

figure(2)
surf(x,y,exact);
title('Exact solution');

figure(3)
plot(log(h), log(error));
title('Order of convergence');
