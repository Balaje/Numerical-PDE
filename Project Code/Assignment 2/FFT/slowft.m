function V = slowft(x)
% Calculating DFT using the slow FFT.
N = size(x,1);

k = 0:1:N-1;
n = k;

% M = exp(-2i*k'*n/N);
M = zeros(N);
for p=1:N
    for q=1:N
        M(p,q) = exp(-1i*2*pi*(k(p))*(n(q))/N);
    end
end

V = M*x;
%disp(t)
end