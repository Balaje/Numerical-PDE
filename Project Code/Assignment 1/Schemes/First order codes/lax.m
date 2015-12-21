% Lax Friedrich Scheme

syms x

x0 = 0;
xf = 1;

N = 100;

h = (xf-x0)/N;
eps = 10^-20;

delt = 0.5*h;

t0 = 0;
tf = 0.5;

x = zeros(N-1,1);
u0 = zeros(N-1,1);

for j=1:N-1
    x(j) = x0+j*h;
    u0(j) = sin(pi*x(j));
end

%plot(x,u0);

M = fix((tf-t0)/delt);

unew = zeros(N-1,1);
for k=1:M
    unew(1) = (u0(2))/2 - 0.5*(delt/h)*(f(u0(2))-f(0)) + eps*(delt/h^2)*(u0(2)-2*u0(1)+0);
    unew(N-1) = u0(N-2)/2 + 0.5*(delt/h)*(f(u0(N-2)))+ eps*(delt/h)*(u0(N-2)-2*u0(N-1));
    for j=2:N-2
        unew(j) = (u0(j+1)+u0(j-1))/2 - 0.5*(delt/h)*(f(u0(j+1))-f(u0(j-1))) + eps*(delt/h^2)*(u0(j+1)-2*u0(j)+u0(j-1));
    end
    u0 = unew;
    plot(x,unew)
    title('Lax Friedrich Scheme');
    pause(0.01)
end

ulax = unew;