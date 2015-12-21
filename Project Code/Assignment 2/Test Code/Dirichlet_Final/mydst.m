%mydst for Dirichlet problem

function V = mydst(U)
[N1,N2]= size(U);

M = zeros(N1);
V = zeros(N1,N2);
for k=1:N1
    for i=1:N1
        M(k,i) = sin(i*k*pi/(N1+1));
    end
end

for j=1:N2
    V(:,j) = M*U(:,j);
end
end