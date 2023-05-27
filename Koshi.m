function [res, countF, trail] = Koshi(f, x0, eps, maxIt)
x = sym('x', [2 1]);
s = matlabFunction(simplifyFraction(transpose(jacobian(sym(f(x)), x))), 'Vars', {x});
deltaCoef = 0.05; trail = {x0};

it = 0; iF = 0; x = x0;
while 1
    trail{end+1} = x;
    sx = s(x); px = x;
    if norm(sx) < eps(1) || it > maxIt(1)
        break;
    end
    [x, i] = GoldSearch(f, -sx, x, norm(sx)*deltaCoef, eps(2), maxIt(2));
    iF = iF + i + 1; it = it + 1;
end
res = x; countF = iF;
end


