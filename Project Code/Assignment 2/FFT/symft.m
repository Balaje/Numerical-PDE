function y = symft(x)
    n = length(x);
    if n == 1
        y = x;
    else
        m = n/2;
        y_top = symft(x(1:2:(n-1)));
        y_bottom = symft(x(2:2:n));
        d = exp(-2 * pi * 1i / n) .^ (0:m-1);
        z = d .* y_bottom;
        y = [ y_top + z , y_top - z ];
    end
end