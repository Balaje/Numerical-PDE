% Program to solve the BVP

% u_t = u_xx + u_yy + f(t,x,y)
% u(x,y,0) = sin(pi*x)*cos(pi*y)
% u(x,0,t) = exp^(-t)*sin(pi*x)
% u(0,y,t) = 0
% u(x,1,t) = -exp(-t)*sin(pi*x)
% u(1,y,t) = 0

clear all
close all
format long

a1 = 0; b1 = 1;
c = 0; d = 1;

t0 = 0; tf = 0.001;

M = 1000; %Time
N = 10; %Space

for q = 1:5
h = (b1-a1)/N;

del_t = (tf-t0)/M;

mu = del_t/(2*h^2);

x = a1:h:b1;
y = c:h:d;

u0 = zeros(N-1,N-1);
unew = u0;

for i=1:N+1
    for j=1:N+1
        u0(i,j) = sin(pi*x(i))*cos(pi*x(j));
    end
end

A = sparse(N-1,N-1);
for j=1:N-1
    A(j,j) = 1+2*mu;
end
for j=2:N-1
    A(j,j-1) = -mu;
    A(j-1,j) = -mu;
end

t = t0;

b = zeros(N-1,1);

%Start Time Loop
for p=1:M
    % Begin x-direction
    for i=1:N+1
        unew(i,1) = exp(-t)*sin(pi*x(i));
        unew(i,N+1) = -exp(-t)*sin(pi*x(i));
    end
    
    for j=2:N %For a fixed j
        for i=2:N
            b(i-1) = mu*u0(i,j-1) + (1-2*mu)*u0(i,j) + mu*u0(i,j+1) + 0.5*del_t*f((t+0.5*del_t),x(i),y(j));
        end
        ut = A\b;
        
        for i=1:N-1
            unew(i+1,j) = ut(i);
        end
    end
    %End x direction
    
    
    % Begin y-direction
    for i=1:N+1
        u0(1,i) = 0;
        u0(N+1,i) = 0;
    end
    
    for i=2:N %For a fixed i
        for j=2:N
            b(j-1) = mu*unew(i,j-1) + (1-2*mu)*unew(i,j) + mu*unew(i,j+1) + 0.5*del_t*f((t+0.5*del_t),x(i),y(j));
        end
        ut = A\b;
        
        for j=1:N-1
            u0(i,j+1) = ut(j);
        end
    end
    %End y direction
    
    t = t+del_t;
end

exact = zeros(N+1,N+1);
for i=1:N+1
    for j=1:N+1
        exact(i,j) = exp(-t)*sin(pi*x(i))*cos(pi*y(j));
    end
end

error(q) = max(max(abs(u0-exact)));
N = N*2;
end

for j=1:q-1
    order(j) = log(error(j)/error(j+1))/log(2);
end
figure(2)
surf(x,y,exact);
figure(1)
surf(x,y,u0);
