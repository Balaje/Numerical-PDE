%mydst for Neumann problem
function V = myidst(U)
[N1,N2] = size(U);

M = zeros(N1);
V = zeros(N1,N2);
for i=1:N1
    for j=1:N1
        M(i,j) = sin(i*(2*j-1)*pi/(2*N1));
    end
end

for j=1:N2
    V(:,j) = (2/N1)*M*U(:,j);
end
end