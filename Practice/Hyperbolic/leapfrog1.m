% Leap frog scheme (Hyperbolic)
% Uses Lax Freidrich

a = 1;
x0 = 0; xf = 2;
N = 100;
h = (xf-x0)/N;

k = 0.96*h;
t0 = 0; tf = 1;
M = fix((tf-t0)/k);

x = x0+h:h:xf-h;
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

mu = k/h;
u = zeros(N-1,1);
u1 = zeros(N-1,1);

% n = 1
u1(1) = u0(2)/2 - a*mu/2*u0(2);
u1(N-1) = u0(N-2)/2 + a*mu/2*u0(N-2);
for j=2:N-2
    u1(j) = (u0(j+1) + u0(j-1))/2 - mu*a/2*(u0(j+1)-u0(j-1));
end

for p=2:M
    u(1) = u0(1) - a*mu*(u1(2));
    for j=2:N-2
        u(j) = u0(j) - a*mu*(u1(j+1)-u1(j-1));
    end
    u(N-1) = u0(N-1) - a*mu*(-u1(N-2));
    
    u0 = u1;
    u1 = u;
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

plot(x,u,x,ue);