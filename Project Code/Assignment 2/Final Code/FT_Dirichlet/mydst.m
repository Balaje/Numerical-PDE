function V = mydst(U)
%mydst for Dirichlet problem

[N1,N2]= size(U); % Get the size of U

M = zeros(N1); % Initialize the 'M' matrix
V = zeros(N1,N2); % Initialize Output 
% Compute M
for k=1:N1
    for i=1:N1
        M(k,i) = sin(i*k*pi/(N1+1));
    end
end
% Obtaining the Fourier Transform
for j=1:N2
    V(:,j) = M*U(:,j);
end
end