% Leap frog scheme (Hyperbolic)
% Uses Godunov scheme

a = -1;
x0 = 0; xf = 2;
N = 100;
h = (xf-x0)/N;

k = 1*h/abs(a);
t0 = 0; tf = 0.25;
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
% u1(1) = u0(2)/2 - a*mu/2*u0(2);
% u1(N-1) = u0(N-2)/2 + a*mu/2*u0(N-2);
% for j=2:N-2
%     u1(j) = (u0(j+1) + u0(j-1))/2 - mu*a/2*(u0(j+1)-u0(j-1));
% end
u1(1) = u0(1) - 0.5*mu*(1+sign(a))*a*(u0(1)) - 0.5*mu*(1-sign(a))*a*(u0(2)-u0(1));
for j=2:N-2
    u1(j) = u0(j) - 0.5*mu*(1+sign(a))*a*(u0(j)-u0(j-1)) - 0.5*mu*(1-sign(a))*a*(u0(j+1)-u0(j));
end
u1(N-1) = u0(N-1) - 0.5*mu*(1+sign(a))*a*(u0(N-1)-u0(N-2)) - 0.5*mu*(1-sign(a))*a*(-u0(N-1));


for p=2:M
    u(1) = u0(1) - a*mu*(u1(2));
    for j=2:N-2
        u(j) = u0(j) - a*mu*(u1(j+1)-u1(j-1));
    end
    u(N-1) = u0(N-1) - a*mu*(-u1(N-2));
    
    u0 = u1;
    u1 = u;
    plot(x,u)
    pause(0.1)
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

plot(x,u);