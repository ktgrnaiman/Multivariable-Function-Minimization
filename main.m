%% Run grid test
f = @(x) (100*(x(2,:) - x(1,:)^2)^2 + (x(1,:) - 1)^2);
xMinA = [1; 2];
%Simulation parameters
eps = [10E-6, 10E-6]; maxIt = [1000, 1000];
%Grid and initial value
gridSize = 10; step = 2; 
x0 = cell(gridSize,gridSize); x0{1,1} = [-5; -5];
%All metrics
itRes = zeros(gridSize,gridSize,3); xRes = cell(gridSize,gridSize,3); fxRes = zeros(gridSize,gridSize,3);
trail = cell(gridSize,gridSize,3); 

for i=1:gridSize
    for j=1:gridSize
        %Filling grid of initial values
        x0{i,j} = x0{1,1} + [step*(i-1);step*(j-1)];
        %Launching algorithms and recording minima, iterations count and trail
        [xRes{i,j,1}, itRes(i,j,1), trail{i,j,1}] = Koshi(f, x0{i,j}, eps, maxIt);
        [xRes{i,j,2}, itRes(i,j,2), trail{i,j,2}] = Newton(f, x0{i,j}, eps(1), maxIt(1));
        [xRes{i,j,3}, itRes(i,j,3), trail{i,j,3}] = Mark(f, x0{i,j}, eps(1), maxIt(1));
        %Finding function value
        for q=1:3
            fxRes(i,j,q) = f(xRes{i,j,q});
        end
    end
end
%% Minima selection
x = sym('x', [2 1]);
jF = matlabFunction(simplifyFraction(transpose(jacobian(sym(f(x)), x))), 'Vars', {x});
h = matlabFunction(simplifyFraction(hessian(sym(f(x)), x)), 'Vars', {x});
rC = 4;

rightResCount = 0; minimaList = [];
for i=1:gridSize
    for j=1:gridSize
        for q=1:3
            exist = 0;
            if round(norm(jF(xRes{i,j,q}))/2, rC) == 0 && det(h(xRes{i,j,q})) > 0
                itRes(i,j,q) = abs(itRes(i,j,q));
                for z=1:length(minimaList)
                    if round(norm(xRes{i,j,q} - minimaList(z).pos), rC) == 0
                        exist1 = 0;
                        for t=1:length(minimaList(z).trail)
                            if norm([i,j] - minimaList(z).trail(t).pos) == 0
                                minimaList(z).trail(t).alg(q) = 1;
                                exist1 = 1;
                            end
                        end
                        if exist1 == 0
                            minimaList(z).trail = [minimaList(z).trail, Trail([i,j], q)];
                        end
                        exist = 1;
                    end
                end
                if exist == 0
                    if isempty(minimaList)
                        minimaList = [Minima(xRes{i,j,q}, fxRes(i,j,q), [i,j], q)];
                    else
                        minimaList = [minimaList, Minima(xRes{i,j,q}, fxRes(i,j,q), [i,j], q)];
                    end
                end
            else
                itRes(i,j,q) = -itRes(i,j,q);
            end
        end
    end
end
minAIndex = 1;
for i=1:length(minimaList)
    if minimaList(i).f < minimaList(minAIndex).f
        minAIndex = i;
    end
end
%% Statistics
right = [0; 0; 0]; wrong = [0; 0; 0];
convSpeed = ones(gridSize,gridSize,3);
genConvSpeed = [1; 1; 1]; genIt = [1; 1; 1];
genElCount = [0; 0; 0]; wrongStepCount = [0; 0; 0];
for i=1:gridSize
    for j=1:gridSize
        for q=1:3
            elCount = 0;
            for p=3:length(trail{i,j,q})-1
                if f(trail{i,j,q}{p}) - f(trail{i,j,q}{p-1}) < 0
                    coef = norm(f(trail{i,j,q}{p}) - f(trail{i,j,q}{p-1})) / norm(f(trail{i,j,q}{p+1}) - f(trail{i,j,q}{p}));
                    if coef ~= Inf && ~isnan(coef) && coef ~= 0
                        convSpeed(i,j,q) = convSpeed(i,j,q) * coef;
                        elCount = elCount + 1;
                    end
                else
                    wrongStepCount(q) = wrongStepCount(q) + 1;
                end
            end
            convSpeed(i,j,q) = convSpeed(i,j,q)^(1/(elCount));
            if itRes(i,j,q) < 0
                wrong(q) = wrong(q) + 1;
            else
                right(q) = right(q) + 1;
                if ~isempty(p)
                    genConvSpeed(q) = genConvSpeed(q) * convSpeed(i,j,q);
                    genIt(q) = genIt(q) * abs(itRes(i,j,q));
                    genElCount(q) = genElCount(q) + 1;
                end
            end
        end
    end
end
for q=1:3
    genConvSpeed(q) = genConvSpeed(q) ^ (1/genElCount(q));
    genIt(q) = genIt(q) ^ (1/genElCount(q));
end
%% Visualization for every starting point
%Plot parameters
plotGridFreq = 51;

for i=1:2:gridSize
    for j=1:2:gridSize
        %Part with values of function preparations. They should be
        %specifically spaced in order to maintain plot representable
        deltaX = x0{i,j}(1) - xRes{i,j,3}(1);
        deltaY = x0{i,j}(2) - xRes{i,j,3}(2);
        X = nonLinspaceMid(xRes{i,j,3}(1) - deltaX, xRes{i,j,3}(1) + deltaX, 1.1, plotGridFreq);
        Y = nonLinspaceMid(xRes{i,j,3}(2) - deltaY, xRes{i,j,3}(2) + deltaY, 1.1, plotGridFreq);
        Z = fGrid(f, X, Y);
        L = nonLinspaceEnd(max(max(Z)), minimaList(minAIndex).f, 1.2, 2*plotGridFreq);
        %Part of contour plot funcionality
        figure('Name', strcat("Contour", num2str(x0{i, j}(1)), num2str(x0{i, j}(2))));
        contour(X, Y, Z, L, 'LineWidth', 0.5); hold on;
        plotTrail(trail{i,j,1},1, 1.2);
        plotTrail(trail{i,j,2},2, 1.2);
        plotTrail(trail{i,j,3},3, 1.2);
        legend("","Koshi","Newton","Markquardt");
        %Part of 3d plot funcionaluty
        figure('Name','3d');
        mesh(X, Y, Z); hold on;
        plotTrail3(f, trail{i,j,1},1, 1.2);
        plotTrail3(f, trail{i,j,2},2, 1.2);
        plotTrail3(f, trail{i,j,3},3, 1.2);
        legend("","Koshi","Newton","Markquardt");
    end
end
hold off;
%% Visualization of minimas
plotGridFreq = 51;
for i = 1:length(minimaList)
    window = figure;   
    hold on;
    X = nonLinspaceMid(x0{1,1}(1), x0{end,1}(1), 1.1, plotGridFreq);
    Y = nonLinspaceMid(x0{1,1}(2), x0{1,end}(2), 1.1, plotGridFreq);
    Z = fGrid(f, X, Y);
    L = nonLinspaceEnd(max(max(Z)), minimaList(minAIndex).f, 1.2, 2*plotGridFreq);
    contour(X, Y, Z, L, 'LineWidth', 0.5);
    for j = 1:length(minimaList(i).trail)
        for z = 1:3
            if minimaList(i).trail(j).alg(z) == 1
                curTrail = trail{minimaList(i).trail(j).pos(1), minimaList(i).trail(j).pos(2), z};
                plot(curTrail{end}(1), curTrail{end}(2), '.','MarkerSize', 20, 'Color', [0, 0, 0]);
                plotTrail(curTrail, z, 1);
            end
        end
    end
    if i == minAIndex
        window.Name = 'Absolute';
    end
    hold off;
end
%% Colormap of expended function calculations
names = ["Koshi", "Newton", "Markquard"];
X = zeros(gridSize,gridSize); Y = zeros(gridSize,gridSize);
for i=1:gridSize
    for j=1:gridSize
        X(i,j) = x0{i,j}(1);
        Y(i,j) = x0{i,j}(2);
    end
end
for i=1:3
    figure('Name',strcat(names(i),"Function calculation count"));
    pcolor(X,Y,abs(itRes(:,:,i))); xlabel("x"); ylabel("y");
    colorbar
    figure('Name',strcat(names(i),"Average convergence speed"));
    pcolor(X,Y,convSpeed(:,:,i)); xlabel("x"); ylabel("y");
    colorbar
end











