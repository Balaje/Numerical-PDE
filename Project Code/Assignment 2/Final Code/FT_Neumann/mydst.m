%mydst for Neumann problem
function V = mydst(U)
[N1, N2] = size(U);

M = zeros(N1);
V = zeros(N1,N2);
for i=1:N1
    for j=1:N1
        M(i,j) = sin(((2*i-1)*j*pi/(2*N1)));
    end
end

M(:,N1) = M(:,N1)/2;
for j=1:N2
    V(:,j) = M*U(:,j);
end

end