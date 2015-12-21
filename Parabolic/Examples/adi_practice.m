% ADI Scheme
% Problem: Ut = Uxx+Uyy + exp(-t)*(2*pi^2-1)*sin(pi*x)*cos(pi*y)
%               U(x,y,0) = sin(pi*x)*cos(pi*y);
%               U = 0 on the veritcal boundary
%               U(x,0) = exp(-t)*sin(pi*x);
%               U(x,1) = -exp(-t)*sin(pi*x);
%

clear all;
close all;

format long

%Assuming same mesh size for y;
x0 = 0; y0 = 0;
x1 = 1; y1 = 1;

t0 = 0;
tf = 1;

N = 10;
M = N;
T = 100;

h = (x1-x0)/N;
del_t = (tf-t0)/M;

mu = del_t/(2*h^2);

x = x0:h:x1;
y = y0:h:y1;
u1 = zeros(N+1);
u2 = zeros(N+1);

for i=1:N+1
    for j=1:N+1
        u1(i,j) = sin(pi*x(i))*cos(pi*y(j));
    end
end

A = zeros(N-1,N-1);

t = t0;
b = zeros(N-1,1);

for k=1:T
    %Begin x direction
    for j=1:N+1
        u2(j,1) = exp(-t)*sin(pi*x(j));
        u2(j,N+1) = -exp(-t)*sin(pi*x(j));
    end
    
    for j=1:N-1
        A(j,j) = 1+2*mu;
    end
    for j=2:N-1
        A(j,j-1) = -mu;
        A(j-1,j) = -mu;
    end
    
    for j=2:M %Fix y
        for i=2:N
            b(i-1) = mu*u1(i,j-1) + (1-2*mu)*u1(i,j) + mu*u1(i,j+1) + del_t*0.5*f((t+0.5*del_t),x(i),y(j));
        end
        ut = A\b;
        for i=1:N-1
            u2(i+1,j) = ut(i);
        end
    end %End x direction
%---------------------------
    
    %Begin y direction
    for j=1:N+1
        u1(1,j) = 0;
        u1(N+1,j) = 0;
    end

    for j=1:N-1
        A(j,j) = 1+2*mu;
    end
    for j=2:N-1
        A(j,j-1) = -mu;
        A(j-1,j) = -mu;
    end
    
    for j=2:N %Fix x
        for i=2:M
            b(i-1) = mu*u2(j-1,i) + (1-2*mu)*u2(j,i) + mu*u2(j+1,i) + del_t*0.5*f((t+0.5*del_t),x(j),y(i));
        end
        ut = A\b;
        for i=1:N-1
            u1(j,i+1) = ut(i);
        end %End y direction
    end
%--------------------------   
%--------------------------
    t = t+del_t;
end

for i=1:N+1
    for j=1:N+1
        exact(i,j) = exp(-t)*sin(pi*x(i))*cos(pi*y(j));
    end
end
figure(1)
surf(x,y,u1);
xlabel('x'); ylabel('y');
figure(2);
surf(x,y,exact);
xlabel('x'); ylabel('y');
