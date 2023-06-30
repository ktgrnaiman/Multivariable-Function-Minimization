function [res, i] = GoldSearch(f, s, x0, h, eps, n)
    ns = norm(s); s = s / ns; g = (1 + sqrt(5))/2;
    eps = sqrt(eps);
    a = 0; b = h; i = 3;

    while f(x0 + s*b) < f(x0 + s*a)
        a = b; b = 2*b; i = i + 1;
    end

    a = a / 2;
    x1 = b - (b - a)/g; x2 = a + (b - a)/g;
    f1 = f(x0 + s*x1); f2 = f(x0 + s*x2);
    while abs(b - a) > eps && i < n
        if f1 >= f2
            a = x1; x1 = x2; x2 = b - (x1 - a); 
            f1 = f2; f2 = f(x0 + s*x2);
        else
            b = x2; x2 = x1; x1 = a + (b - x2); 
            f2 = f1; f1 = f(x0 + s*x1);
        end
        i = i + 1;
    end
    res = x0 + s*((a + b)/2);
end

