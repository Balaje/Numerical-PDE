% Einguist Osher Scheme
% ut + (u^2/2)x = eps*uxx;

clear all
x0 = 0;
xf = 1;
t0 = 0;
tf = 1;
eps = 0.09;

[x1,exact] = genexact(x0,xf,t0,tf,eps);

N = 10;
z = 4;
order = zeros(z-1,1);
error2 = zeros(z,1);
h = zeros(z,1);

for p=1:z

h(p) = (xf-x0)/N;
delt = h(p)/(1+2*eps/h(p));

lambda = delt/h(p);
mu = delt/(h(p)^2);

x = zeros(N-1,1);
u0 = zeros(N-1,1);

for j=1:N-1
    x(j) = x0+j*h(p);
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
    %pause(0.01);
    u0 = unew;
    t = t+delt;
end

    error1 = zeros(gcd(size(x1,1)+1,size(x,1)+1),1);
    for i=1:N-1
        for j=1:size(exact,1)-1
            if(x(i) == x1(j))
                error1(i) = abs(exact(j)-unew(i));
            end
        end
    end
    error2(p) = max(error1);
    N = N*2;
end
for j=1:p-1
    order(j) = log(error2(j)/error2(j+1))/log(2);
end
