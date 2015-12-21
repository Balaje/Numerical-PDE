% Central Difference
% ut + (u^2/2)x = eps*uxx;

syms x

x0 = 0;
xf = 5;

eps = 0.01;

h = 0.01;
delt = h^2/(4*eps);
mu = delt/h;
tau = delt/h^2;

t0 = 0;
tf = 5;

N = fix((xf-x0)/h);
x = zeros(N-1,1);
u0 = zeros(N-1,1);

for j=1:N-1
    x(j) = x0+j*h;
    if(x(j) <= 1)
        u0(j) = 3*sin(pi*x(j));
    end
    if(x(j) > 1)
        u0(j) = 0;
    end
end

%plot(x,u0);

M = fix((tf-t0)/delt);

unew = zeros(N-1,1);
for k=1:M
    unew(1) = u0(1) - 0.5*mu*(f(u0(2))-f(0)) + eps*tau*(u0(2)-2*u0(1)+0);
    unew(N-1) = u0(N-1) + 0.5*mu*(f(u0(N-2))) + eps*tau*(u0(N-2)-2*u0(N-1));
    for j=2:N-2
        unew(j) = u0(j) - 0.5*mu*(f(u0(j+1))-f(u0(j-1))) + eps*tau*(u0(j+1)-2*u0(j)+u0(j-1));
    end
    u0 = unew;
    plot(x,unew)
    axis([x0,xf,0,3]);
    title('Central Difference Scheme');
    pause(0.0001)
end

%plot(x,unew)
unl = unew;