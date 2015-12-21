% Godunov Scheme
% ut + (u^2/2)x = eps*uxx;

x0 = 0;
xf = 5;

N = 100;

h = (xf-x0)/N;
eps = 0.01;

delt = 0.5*h;

lambda = delt/h;
mu = delt/(h^2);

t0 = 0;
tf = 5;

x = zeros(N-1,1);
u0 = zeros(N-1,1);
for j=1:N-1
    x(j) = x0+j*h;
    if(x(j) < 1)
        u0(j) = sin(pi*x(j));
    end
    if(x(j) > 1)
        u0(j) = 0;
    end
end
u0g = u0;

M = fix((tf-t0)/delt);

unew = zeros(N-1,1);
t = t0;
for k=1:M
    unew(1) = u0(1) - lambda*(fluxfunction(u0(1),u0(2)) - fluxfunction(0,u0(1))) + mu*eps*(u0(2) - 2*u0(1) + 0);
    unew(N-1) = u0(N-1) - lambda*(fluxfunction(u0(N-1),0) - fluxfunction(u0(N-2),u0(N-1))) + mu*eps*(u0(N-2) - 2*u0(N-1) + 0);
    for j=2:N-2
        unew(j) = u0(j) - lambda*(fluxfunction(u0(j),u0(j+1)) - fluxfunction(u0(j-1),u0(j))) + mu*eps*(u0(j+1) - 2*u0(j) + u0(j-1));
    end
    plot(x,unew)
    title('Godunov Scheme');
    ylim([0,1]);
    pause(0.01);
    u0 = unew;
    t = t+delt;
end

ugod = unew;