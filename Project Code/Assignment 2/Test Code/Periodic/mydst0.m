function V = mydst0(U)
[N1, N2] = size(U);

M = zeros(N1);
V = zeros(N1,N2);

for k=1:N1
    for i=1:N1
        M(k,i) = (2/(N1))*cos(2*(k-1)*(i-1)*pi/(N1));
    end
end


for j=1:N2
    V(:,j) = M*U(:,j);
end

end