f = @(x) (1 + 1/4000*(x(1,:)^2 + x(2,:)^2) - cos(x(1,:))*cos(x(2,:)/sqrt(2)));
eps = [10E-6, 10E-6]; maxIt = [1000, 1000]; x0 = [-10; -10];
it = [0 0 0]; x = [0 0 0];
hold on;

[x1, it(1)] = Koshi(f, x0, eps, maxIt);
[x2, it(2)] = Newton(f, x0, eps(1), maxIt(1));
[x3, it(3)] = Mark(f, x0, eps(1), maxIt(1));
res = [f(x1), f(x2), f(x3)];

delta = norm(x0 - x);
X = linspace(x3(1) - delta, x3(1) + delta, 200);
Y = linspace(x3(2) - delta, x3(2) + delta, 200);
Z = fGrid(f, X, Y);
contour(X, Y, Z, logspace(-100,100,1000));
str = {'\color{blue} Cauchy', '\color{red} Newton', '\color{green} Marquard'};
legend(str, Interpreter ='tex');

