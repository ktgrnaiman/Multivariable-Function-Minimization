%%Run grid test
f = @(x) (4*x(1,:)^2 - 2.1*x(1,:)^4 + 1/3*x(1,:)^6 + x(1,:)*x(2,:) - 4*x(2,:)^2 + 4*x(2,:)^4);
eps = [10E-6, 10E-6]; maxIt = [1000, 1000]; ifPrint = 0;
it = zeros(10,10,3); xMin = cell(10,10,3); fxMin = zeros(10,10,3);
x0 = cell(10,10); x0{1,1} = [-10; -10]; step = 1; stepCount = 20;
hold on;

for i=1:stepCount
    for j=1:stepCount
        if i==3 && j==3
            ifPrint = 1;
        end
        x0{i,j} = x0{1,1} + [step*(i-1);step*(j-1)];
        [xMin{i,j,1}, it(i,j,1)] = Koshi(f, x0{i,j}, eps, maxIt, ifPrint);
        [xMin{i,j,2}, it(i,j,2)] = Newton(f, x0{i,j}, eps(1), maxIt(1), ifPrint);
        [xMin{i,j,3}, it(i,j,3)] = Mark(f, x0{i,j}, eps(1), maxIt(1), ifPrint);
        for q=1:3
            fxMin(i,j,q) = f(xMin{i,j,q});
        end
        if ifPrint == 1
            ifPrint = 0;
            delta = 5*norm(x0{i,j} - xMin{i,j,3});
            X = nonLinspaceMid(xMin{i,j,3}(1) - delta, xMin{i,j,3}(1) + delta, 1.1, 50);
            Y = nonLinspaceMid(xMin{i,j,3}(2) - delta, xMin{i,j,3}(2) + delta, 1.1, 50);
            L = nonLinspaceMid(-1000,1000, 1.1, 100);
            Z = fGrid(f, X, Y);
            contour(X, Y, Z, L, 'LineWidth', 0.5);
            str = {'\color{blue} Cauchy', '\color{red} Newton', '\color{green} Marquard'};
            legend(str, Interpreter ='tex');
        end
    end
end

hold off;
%% Visualize
rightResCount = 0;
X = zeros(10,10); Y = zeros(10,10);
for i=1:stepCount
    for j=1:stepCount
        X(i,j) = x0{i,j}(1);
        Y(i,j) = x0{i,j}(2);
    end
end
for i=1:3
    figure
    s = pcolor(X,Y,it(:,:,i)); xlabel("x"); ylabel("y");
    colorbar
end
for i=1:10
    for j=1:10
        for q=1:3
            if round(xMin{i,j,q}(1),3) == round(0.089842,3) && round(xMin{i,j,q}(2),3) == round(-0.712656,3)
                rightResCount = rightResCount + 1;
            end
        end
    end
end






