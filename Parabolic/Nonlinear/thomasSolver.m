function V = thomasSolver(A,F)

N = size(F,1);
a = zeros(N,1);
b = zeros(N,1);
c = zeros(N,1);

a(1) = 0;
for j=1:N
    b(j) = A(j,j);
end
for j=2:N
    a(j) = A(j,j-1);
    c(j-1) = A(j-1,j);
end
c(N) = 0;

% Thomas Algorithm

for k=2:N
    m = a(k)/b(k-1);
    b(k) = b(k) - m*c(k-1);
    F(k) = F(k) - m*F(k-1); 
end

V = zeros(N,1);

V(N) = F(N)/b(N);

for k=N-1:-1:1
    V(k) = (F(k) - c(k)*V(k+1))/b(k);
end

end