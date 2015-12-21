% Godunov scheme for hyperbolic PDEs
clear
close
a = 1;
x0 = 0; xf = 2;
N = 200;

for q=1:5
h = (xf-x0)/N; 

t0 = 0; tf = 0.1;
k = 0.5*h/abs(a);
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
    u(1) = u0(1) - 0.5*mu*(1+sign(a))*a*(u0(1)) - 0.5*mu*(1-sign(a))*a*(u0(2)-u0(1));
    for j=2:N-2
        u(j) = u0(j) - 0.5*mu*(1+sign(a))*a*(u0(j)-u0(j-1)) - 0.5*mu*(1-sign(a))*a*(u0(j+1)-u0(j));
    end
    u(N-1) = u0(N-1) - 0.5*mu*(1+sign(a))*a*(u0(N-1)-u0(N-2)) - 0.5*mu*(1-sign(a))*a*(-u0(N-1));
   
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
    elseif(x(j)-a*tf >= 1)
        ue(j) = 0;
    end
end

error(q) = max(abs(u-ue));

N = N*2;
end
%plot(x,u)

for j=1:q-1
    order(j) = log(error(j)/error(j+1))/log(2);
end
plot(x,u,x,ue)