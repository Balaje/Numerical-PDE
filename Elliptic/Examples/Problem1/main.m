% Problem 3 Lab Sheet - Elliptic 
clear all
close all
format long

x0 = 0; xM = 2;
y0 = 0; yN = 1;

M = 40;
N = 80;

hx = (xM-x0)/M;
hy = (yN-y0)/N;

hx1 = hx^2;
hy1 = hy^2;

x = zeros(M-1,1);
for i=1:M-1
    x(i) = x0+i*hx;
end

y = zeros(N-1,1);
for i=1:N-1
    y(i) = x0+i*hy;
end

M1 = (M-1)*(N-1);
A = zeros(M1);
F = zeros(M1,1);

F(1) = f(x(1),y(1)) - (1/hx1 + 1/hy1);
for i=2:M-2
    k = i;
    F(k) = f(x(i),y(1)) - 1/hy1;
end
F(M-1) = f(x(M-1),y(1)) - (exp(2*y(1))/hx1 + 1/hy1);

for j=2:N-2
    F(1+(M-1)*(j-1)) = f(x(1),y(j)) - 1/hx1;
    for i=2:M-2
        k = i+(M-1)*(j-1);
        F(k) = f(x(i),y(j));
    end
    F(M-1 + (M-1)*(j-1)) = f(x(M-1),y(j)) - (exp(2*y(j))/hx1);
end

F(1+(M-1)*(N-2)) = f(x(1),y(N-1)) - (1/hx1 + exp(x(1))/hy1);
for i=2:M-2
    k = i+(M-1)*(N-2);
    F(k) = f(x(i),y(N-1)) - exp(x(i))/hy1;
end
F(M-1+(M-1)*(N-2)) = f(x(M-1),y(N-1)) - (exp(2*y(N-1))/hx1 + exp(x(M-1))/hy1);


for j = 1:N-1
  for i=1:M-1
    k = i + (j-1)*(M-1);
    A(k,k) = -2/hx1 -2/hy1;
    if i == 1
        A(k,k+1) = 1/hx1;
    else
       if i==M-1
         A(k,k-1) = 1/hx1;
       else
        A(k,k-1) = 1/hx1;
        A(k,k+1) = 1/hx1;
       end
    end

%-- y direction --------------

    if j == 1
        A(k,k+M-1) = 1/hy1;
    else
       if j==N-1
         A(k,k-(M-1)) = 1/hy1;
       else
          A(k,k-(M-1)) = 1/hy1;
          A(k,k+M-1) = 1/hy1;
       end
     end

  end
end

U = A\F;

j = 1;
for k=1:M1
  i = k - (j-1)*(M-1) ;
  Uapp(i,j) = U(k);
  Uex(i,j) = exp(x(i)*y(j));
  j = fix(k/(M-1)) + 1;
end

e = max( max( abs(Uapp-Uex)))       % The maximum error

surf(y,x,Uapp); 
title('The Approximate solution plot'); 
xlabel('x'); 
ylabel('y');

figure(2); 
surf(y,x,Uex); 
title('The Exact solution plot'); 
xlabel('x');
ylabel('y');
