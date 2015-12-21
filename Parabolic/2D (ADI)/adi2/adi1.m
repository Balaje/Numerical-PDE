% ADI.

a = 0; b = 1;
c = 0; d = 1;
tf = 1; t0 = 0;

N = 100;
M = 100;

h = (b-a)/N;
k = (tf-t0)/M;

mu = k/(2*h^2);

x = a+h:h:b-h;
y = c+h:h:d-h;

U = zeros(N-1);
U12 = U;
U1 = U12;

for i=1:N-1
    for j=1:N-1
        U(i,j) = sin(pi*x(i))*cos(pi*y(j));
    end
end

A = zeros(N-1);
A(1,1) = 1+2*mu;
for j=2:N-1
    A(j,j) = 1+2*mu;
    A(j-1,j) = -mu;
    A(j,j-1) = -mu;
end
F = zeros(N-1,1);
t = 0;
for p=1:M
    %j=1
    for i=1:N-1
        F(i) = U(i,1) + mu*(exp(-t)*sin(pi*x(i)) - 2*U(i,1) + U(i,2) ) +...
            k*0.5*f1(x(i),y(1),t+0.5*k);
    end
    U12(:,1) = A\F;
    % j = 2:N-2
    for j=2:N-2
        for i=1:N-1
            F(i) = U(i,j) + mu*(U(i,j-1) - 2*U(i,j) + U(i,j+1) ) +...
                k*0.5*f1(x(i),y(j),t+0.5*k);
        end
        U12(:,j) = A\F;
    end
    %j=N-1
    for i=1:N-1
        F(i) = U(i,N-1) + mu*(-exp(-t)*sin(pi*x(i)) - 2*U(i,N-1) + U(i,N-2) ) +...
            k*0.5*f1(x(i),y(N-1),t+0.5*k);
    end
    %-------------------------------------------------
     %i=1
     F(1) = U12(1,1) + mu*(exp(-(t+0.5*k))*sin(pi*x(1)) - 2*U12(1,1) + U12(2,1) ) +...
            k*0.5*f1(x(1),y(1),t+0.5*k);
     F(N-1) = U12(1,N-1) + mu*(exp(-(t+0.5*k))*sin(pi*x(N-1)) - 2*U12(1,N-1) + U12(2,N-1) ) +...
            k*0.5*f1(x(N-1),y(N-1),t+0.5*k);   
    for j=2:N-2
        F(j) = U12(1,j) + mu*(exp(-(t+0.5*k))*sin(pi*x(1)) - 2*U12(1,j) + U12(2,j) ) +...
            k*0.5*f1(x(1),y(j),t+0.5*k);
    end
    U1(:,1) = A\F;
    % i = 2:N-2
    for i=2:N-2
        for j=1:N-1
            F(j) = U12(i,j) + mu*(U12(i-1,j) - 2*U12(i,j) + U12(i+1,j) ) +...
                k*0.5*f1(x(i),y(j),t+0.5*k);
        end
        U1(:,j) = A\F;
    end
    %i=N-1
    F(1) = U12(N-1,1) + mu*(exp(-(t+0.5*k))*sin(pi*x(1)) - 2*U12(N-1,1) + U12(N-2,1) ) +...
            k*0.5*f1(x(N-1),y(1),t+0.5*k);
    F(N-1) = U12(N-1,N-1) + mu*(exp(-(t+0.5*k))*sin(pi*x(N-1)) - 2*U12(N-1,N-1) + U12(N-2,N-1) ) +...
            k*0.5*f1(x(N-1),y(N-1),t+0.5*k); 
    for j=2:N-2
        F(j) = U12(N-1,j) + mu*(-exp(-(t+0.5*k))*sin(pi*x(N-1)) - 2*U12(N-1,j) + U12(N-2,j) ) +...
            k*0.5*f1(x(N-1),y(j),t+0.5*k);
    end
    t = t+k;
    
    U = U1;
end

mesh(x,y,U1)