function V = f(t,x)
    V = exp(-t)*sin(pi*x)*(pi^2-exp(-t)*sin(pi*x));
end

