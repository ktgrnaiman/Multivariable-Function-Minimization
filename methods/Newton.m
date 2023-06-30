function [res, cF, trail] = Newton(f, x0, eps, maxIt)
x = sym('x', [2 1]);
j = matlabFunction(simplifyFraction(transpose(jacobian(sym(f(x)), x))), 'Vars', {x});
h = matlabFunction(simplifyFraction(hessian(sym(f(x)), x)), 'Vars', {x});

x = x0; i = 0; iF = 0;
trail = {x0};

while 1
    trail{end+1} = x;
    jx = j(x); 
    if norm(jx) < eps || i > maxIt
        break;
    end
    hx = h(x);
    px = x; x = x - hx\jx;
    i = i + 1; iF = iF + 2;
end
res = x; cF = iF;
end


