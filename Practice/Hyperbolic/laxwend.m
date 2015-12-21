% Lax wendroff Lab exam T2

N = 200;
x0 = 0; xf = 2;
h = (xf-x0)/N;

t0 = 0; tf = 1;

k = 0.8*(h);

M = fix((tf-t0)/k);

x = x0+h:h:xf-h;

lambda = k/h;

u0 = zeros(N-1,1);
for j=1:N-1
    if(x(j) < 0)
        u0(j) = 0;
    elseif(0 <= x(j) && x(j) < 0.5)
        u0(j) = x(j);
    elseif(0.5 <= x(j) && x(j) < 1)
        u0(j) = 1-x(j);
    elseif(x(j) >= 1)
        u0(j) = 0;
    end
end

u = zeros(N-1,1);

for p=1:M
    u(1) = u0(1) - 0.5*lambda*(u0(2)) + 0.5*lambda^2*(u0(2)-2*u0(1));
    for j=2:N-2
        u(j) = u0(j) - 0.5*lambda*(u0(j+1)-u0(j-1)) + 0.5*lambda^2*(u0(j+1)-2*u0(j)+u0(j-1));
    end
    u(N-1) = u0(N-1) - 0.5*lambda*(-u0(N-2)) + 0.5*lambda^2*(u0(N-2)-2*u0(N-1));
    
    u0 = u;
end

exact = zeros(N-1,1);
for j=1:N-1
    if(x(j)-tf < 0)
        exact(j) = 0;
    elseif(0 <= x(j)-tf && x(j)-tf < 0.5)
        exact(j) = x(j)-tf;
    elseif(0.5 <= x(j)-tf && x(j)-tf < 1)
        exact(j) = 1-x(j)+tf;
    elseif(x(j)-tf >= 1)
        exact(j) = 0;
    end
end

plot(x,exact,x,u);
