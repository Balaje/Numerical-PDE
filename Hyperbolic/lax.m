% Advection equation LAx Friedrich

t0 = 0;
tf = 1;
x0 = -1;
xf = 5;
a = 1;

N = 300;

h = (xf-x0)/300;

delt = 0.5*h;

x = zeros(N-1,1);
u0 = zeros(N-1,1);
for j=1:N-1
    x(j) = x0+j*h;
    if(x(j) < 0)
        u0(j) = 1;
    end
    if(x(j) >=0 && x(j) <=1)
        u0(j) = 2*x(j)^3 -3*x(j)^2+1;
    end
    if(x(j)>1)
        u0(j) = 0;
    end
end

M = fix((tf-t0)/delt);

unew = zeros(N-1,1);
for k=1:M
    unew(N-1) = (u0(N-2))/2 + (a*delt/(2*h))*(u0(N-2));
    unew(1) = (u0(2)+1)/2 - (a*delt/(2*h))*(u0(2)-1);
    for j=2:N-2
        unew(j) = (u0(j+1)+u0(j-1))/2 - (a*delt/(2*h))*(u0(j+1)-u0(j-1));
    end
    plot(x,unew);
    pause(0.01);
    u0 = unew;
end

exact = zeros(N-1,1);
for j=1:N-1
    if(x(j) < tf)
        exact(j) = 1;
    end
    if(x(j) >= tf && x(j) <= 1+tf)
        exact(j) = 2*(x(j)-tf)^3 - 3*(x(j)-tf)^2 + 1;
    end
    if(x(j)> 1+tf)
        exact(j) = 0;
    end
end
plot(x,unew,x,exact);