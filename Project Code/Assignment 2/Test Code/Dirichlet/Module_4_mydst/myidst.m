function V = myidst(U)
N = size(U,1);

M = zeros(N);
V = zeros(N);
for i=1:N
    for j=1:N
        M(i,j) = sin(i*j*pi/(N+1));
    end
end

for j=1:N
    V(:,j) = 2/(N+1)*M*U(:,j);
end
end