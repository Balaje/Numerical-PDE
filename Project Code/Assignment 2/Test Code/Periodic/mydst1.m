function V = mydst1(U)
[N1, N2] = size(U);

V = zeros(N1,N2);
%{
for j=2:N2
    V(:,j) = V(:,j) + (2/N1)*sin(2*pi/N1)*U(:,j);
end
%}

M = zeros(N1,N1-1);
for i=1:N1
    for j=1:N1-1
        M(i,j) = (2/N1)*sin(2*(i-1)*j*pi/N1);
    end
end

for j=1:N2
    V(:,j) = M*U(2:N1,j);
end
end