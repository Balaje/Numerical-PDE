function [x,exact] = genexact(x0,xf,t0,tf,eps)

N = 700;

h = (xf-x0)/N;
delt = h/(1+2*eps/h);

lambda = delt/h;
mu = delt/(h^2);

x = zeros(N-1,1);
u0 = zeros(N-1,1);

for j=1:N-1
    x(j) = x0+j*h;
    if(x(j) <= 1)
        u0(j) = sin(pi*x(j));
    end
    if(x(j) > 1)
        u0(j) = 0;
    end
end

%plot(x,u0);

M = fix((tf-t0)/delt);

unew = zeros(N-1,1);
t = t0;
for k=1:M+1
    unew(1) = u0(1) - lambda*(Flux(u0(1),u0(2)) - Flux(0,u0(1))) + mu*eps*(u0(2) - 2*u0(1) + 0);
    unew(N-1) = u0(N-1) - lambda*(Flux(u0(N-1),0) - Flux(u0(N-2),u0(N-1))) + mu*eps*(u0(N-2) - 2*u0(N-1) + 0);
    for j=2:N-2
        unew(j) = u0(j) - lambda*(Flux(u0(j),u0(j+1)) - Flux(u0(j-1),u0(j))) + mu*eps*(u0(j+1) - 2*u0(j) + u0(j-1));
    end
    %plot(x,unew)
    %axis([x0,xf,0,1]);
    %title('Einguist Osher Scheme');
    %spause(0.01);
    u0 = unew;
    t = t+delt;
end

exact = unew;
end

