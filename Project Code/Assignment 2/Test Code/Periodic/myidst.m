function V = myidst(U0,U1)
[N1,N2] = size(U0);

M0 = zeros(N1,N1/2+1);
M1 = M0;
V = zeros(N1,N2);
for i=1:N1
    for j=1:N1/2+1
        M0(i,j) = cos(2*(i-1)*(j-1)*pi/(N1));
        M1(i,j) = sin(2*(i-1)*(j-1)*pi/(N1));
    end
end

M0(:,1) = M0(:,1)/2;
M0(:,N1/2+1) = M0(:,N1/2+1)/2;
M1(:,1) = M1(:,1)/2;
M1(:,N1/2+1) = M1(:,N1/2+1)/2;

for j=1:N2
    V(:,j) = M0*U0(1:N1/2+1,j)+M1*U1(1:N1/2+1,j);
end
end