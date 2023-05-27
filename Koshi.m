function [res, countF] = Koshi(f, x0, eps, maxIt, ifPrint)
x = sym('x', [2 1]);
s = matlabFunction(simplifyFraction(transpose(jacobian(sym(f(x)), x))), 'Vars', {x});
deltaCoef = 0.05;

it = 0; iF = 0; x = x0;
while 1
    sx = s(x); px = x;
    if norm(sx) < eps(1) || it > maxIt(1)
        break;
    end
    [x, i] = GoldSearch(f, -sx, x, norm(sx)*deltaCoef, eps(2), maxIt(2));
    iF = iF + i + 1; it = it + 1;
end
res = x; countF = iF;
end


