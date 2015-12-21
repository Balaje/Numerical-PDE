% An unstable scheme (FTCS)

N = 100;
M = 100;
b = 2; a = 0;

h = (b-a)/N; 

t0 = 0; tf = 5;
k = (tf-t0)/M;
u0 = zeros(N-1,1);
x = a+h:h:b-h;

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

mu = k/(2*h);
u = zeros(N-1,1);

for p=1:M
    u(1) = u0(1) - mu*(u0(2));
    u(N-1) = u0(N-1) - mu*(-u0(N-2));
    for j=2:N-2
        u(j) = u0(j) - mu*(u0(j+1)-u0(j-1));
    end
    
    u0 = u;
end
ue = zeros(N-1,1);
for j=1:N-1
    if(x(j)-tf < 0)
        ue(j) = 0;
    elseif(0 <= x(j)-tf && x(j)-tf < 0.5)
        ue(j) = x(j)-tf;
    elseif(0.5 <= x(j)-tf && x(j)-tf < 1)
        ue(j) = 1-x(j)+tf;
    elseif(x(j)-tf >= 1)
        ue(j) = 0;
    end
end

plot(x,u)