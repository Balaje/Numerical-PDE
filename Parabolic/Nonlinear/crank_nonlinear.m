% Crank Nicolson for non linear problem

clear all;
close all;

format long;

a = 0;
b = 1;
t0 = 0;
tf = 1;

N = 10;
M = 1000;

del_t = (tf-t0)/M;

for p=1:5
h(p) = (b-a)/N;

mu = del_t/(2*h(p)^2);

u0 = zeros(N-1,1);
x = zeros(N-1,1);
for j=1:N-1
    x(j) = j*h(p);
    u0(j) = sin(pi*x(j));
end

F = zeros(N-1,1);
A = zeros(N-1,N-1);

U1 = zeros(N-1,1);
U = ones(N-1,1);
error = 1000;
tol = 10^-8;

t=t0;
while error > tol
    for k=1:M
        F(1) = (1+2*mu-del_t*0.5*u0(1)+del_t*0.5)*U(1) - del_t*0.25*(U(1))^2 ...
            - mu*U(2) + (2*mu-1)*u0(1) - mu*u0(2) - del_t*(u0(1))^2*0.25 + ...
            del_t*0.5*u0(1) - del_t*f((t+0.5*del_t),x(1));
        for j=2:N-2
            F(j) = -mu*U(j-1) + (1+2*mu-0.5*del_t*u0(j)+0.5*del_t)*U(j) - 0.25*del_t*(U(j))^2 ...
                - mu*U(j+1) - mu*u0(j-1) + (2*mu-1)*u0(j) - mu*u0(j+1) - ...
                del_t*(u0(j))^2*0.25 + del_t*0.5*u0(j) - del_t*f((t+0.5*del_t),x(j));
        end
         F(N-1) = -mu*U(N-2) + (1+2*mu-del_t*0.5*u0(1)+del_t*0.5)*U(N-1) - del_t*0.25*(U(N-1))^2 ...
            + (2*mu-1)*u0(N-1) - mu*u0(N-2) - del_t*(u0(N-1))^2*0.25 + ...
            del_t*0.5*u0(N-1) - del_t*f((t+0.5*del_t),x(N-1));
        
        % Jacobian
        for j=1:N-1
            A(j,j) = 1+2*mu-del_t*0.5*u0(j)+del_t*0.5-del_t*0.5*U(j);
        end
        for j=2:N-1
            A(j,j-1) = -mu;
            A(j-1,j) = -mu;
        end
        DELTA = A\F;
        U1 = U - DELTA;
        error = max(abs(U1-U));
        U = U1;
        u0 = U;
        t = t+del_t;
    end
end

exact = zeros(N-1,1);
for j=1:N-1
    exact(j) = exp(-t)*sin(pi*x(j));
end

error1(p) = max(abs(U1-exact));

N = N*2;
end
for k=1:p-1
    order(k) = log(error1(k)/error1(k+1))/log(2);
end

figure(1);
plot(x,exact,x,U1,'-*');

figure(2);
plot(log(h),log(error1));