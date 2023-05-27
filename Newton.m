function [res, cF] = Newton(f, x0, eps, maxIt, ifPrint)
x = sym('x', [2 1]);
j = matlabFunction(simplifyFraction(transpose(jacobian(sym(f(x)), x))), 'Vars', {x});
h = matlabFunction(simplifyFraction(hessian(sym(f(x)), x)), 'Vars', {x});

hold on;
x = x0; i = 0; iF = 0;
while 1
    jx = j(x); 
    if norm(jx) < eps || i > maxIt
        break;
    end
    hx = h(x);
    px = x; x = x - hx\jx;
    i = i + 1; iF = iF + 2;
    if ifPrint == 1
        plot([px(1), x(1)], [px(2), x(2)], 'r', 'LineWidth', 1);
    end
end
res = x; cF = iF;
end


