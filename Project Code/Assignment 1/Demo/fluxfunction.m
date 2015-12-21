function V = fluxfunction(a,b)
% Godunov
V = max(f(max(a,0)),f(min(b,0)));

end

