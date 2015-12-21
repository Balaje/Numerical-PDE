% Backward moving wave

a = -1;
x0 = 0; xf = 2;
N = 200;
h = (xf-x0)/N; 

t0 = 0; tf = 0.1;
k = 0.9*h/abs(a);
M = fix((tf-t0)/k);

u0 = zeros(N-1,1);
x = x0+h:h:xf-h;

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

mu = k/h;
u = zeros(N-1,1);
for p=1:M
    for j=1:N-2
        u(j) = u0(j) - a*mu*(u0(j+1)-u0(j));
    end
    u(N-1) = u0(N-1) - a*mu*(-u0(N-1));
    
    u0 = u;
end
ue = zeros(N-1,1);
for j=1:N-1
    if(x(j)-a*tf < 0)
        ue(j) = 0;
    elseif(0 <= x(j)-a*tf && x(j)-a*tf < 0.5)
        ue(j) = x(j)-a*tf;
    elseif(0.5 <= x(j)-a*tf && x(j)-a*tf < 1)
        ue(j) = 1-x(j)+a*tf;
    elseif(x(j)-a*tf >= 0)
        ue(j) = 0;
    end
end

plot(x,u,x,ue)