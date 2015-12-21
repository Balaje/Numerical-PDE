%myidst for Dirichlet problem

function V = myidst(U)
% U is a N1 x N2 array
[N1,N2] = size(U);

M = zeros(N1); % Initialize the M matrix
V = zeros(N1,N2); % Initialize Output
for i=1:N1
    for k=1:N1
        M(i,k) = sin(i*k*pi/(N1+1)); % Compute M
    end
end

% Compute Inverse Fourier Transform
for j=1:N2
    V(:,j) = 2/(N1+1)*M*U(:,j);
end
end