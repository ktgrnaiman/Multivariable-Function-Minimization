function [res, cF, trail] = Mark(f, x0, eps, maxIt)
x = sym('x', [2 1]);
j = matlabFunction(simplifyFraction(transpose(jacobian(sym(f(x)), x))), 'Vars', {x});
h = matlabFunction(simplifyFraction(hessian(sym(f(x)), x)), 'Vars', {x});

it = 0; iF = 2; x = x0; trail = {x0};
l = 10E4; jx = j(x); hx = h(x);

while norm(jx) > eps && it < maxIt
    trail{end+1} = x;
    s = - inv(hx + l*eye(size(jx, 1)))*jx;
    px = x; x = x + s;
    if f(x) < f(px)
        l = l/2; jx = j(x); hx = h(x);
        iF = iF + 3;
    else
        x = px; l = 2*l; iF = iF + 3;
    end
    it = it + 1;
end
res = x; cF = iF;
end

