%myidst for Dirichlet problem

function V = myidst(U)
[N1,N2] = size(U);

M = zeros(N1);
V = zeros(N1,N2);
for i=1:N1
    for k=1:N1
        M(i,k) = sin(i*k*pi/(N1+1));
    end
end

for j=1:N2
    V(:,j) = 2/(N1+1)*M*U(:,j);
end
end