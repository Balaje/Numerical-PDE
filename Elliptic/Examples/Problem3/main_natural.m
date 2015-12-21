% Program to solve problem 5 using natural ordering

a = 0;
b = 1;
c = 0;
d = 1;

N = 4;
M = 4;

hx = (b-a)/N;
hx1 = hx*hx;

hy = (d-c)/M;
hy1 = hy*hy;

x = zeros(N-1,1);
y = zeros(M+1,1);

M1 = (N-1)*(M+1);

A = zeros(M1);
F = zeros(M1,1);

%k=i+(j-1)*(N-1)
%j=1
F(1) = f(x(1),y(1)) - (0.5*cos(2*pi*y(1)))/hx1;
for i=2:N-2
    k=i;
    F(k) = f(x(i),y(1));
end
F(N-1) = f(x(N-1),y(1)) - 1/hx1;

for j=2:M
    F(1+(j-1)*(N-1)) = f(x(1),y(j)) - (0.5*cos(2*pi*y(j)))/hx1;
    for i=2:N-2
        k = i+(j-1)*(N-1);
        F(k) = f(x(i),y(j));
    end
    F(N-1+(j-1)*(N-1)) = f(x(N-1),y(j)) - 1/hx1;
end

%j=M+1
F(1+M*(N-1)) = f(x(1),y(M+1)) - (0.5*cos(2*pi*y(M+1)))/hx1;
for i=2:N-2
    k = i+M*(N-1);
    F(k) = f(x(i),y(j));
end
F(N-1+M*(N-1)) = f(x(N-1),y(M+1)) - 1/hx1;

for j = 1:M
  for i=1:N-1
    k = i + (j-1)*(N-1);
    A(k,k) = -2/hx1 -2/hy1;
    if i == 1
        A(k,k+1) = 1/hx1;
    else
       if i==N-1
         A(k,k-1) = 1/hx1;
       else
        A(k,k-1) = 1/hx1;
        A(k,k+1) = 1/hx1;
       end
    end

%-- y direction --------------

    if j == 1
        A(k,k+N-1) = 2/hy1;
    else
       if j==M
         A(k,k-(N-1)) = 2/hy1;
       else
          A(k,k-(N-1)) = 1/hy1;
          A(k,k+N-1) = 1/hy1;
       end
     end

  end
end

U = A\F;

j = 1;
for k=1:M1
  i = k - (j-1)*(M-1) ;
  Uapp(i,j) = U(k);
  %Uex(i,j) = exp(x(i)*y(j));
  j = fix(k/(M-1)) + 1;
end