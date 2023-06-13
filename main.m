%% Run grid test
f = @(x) (4*x(1,:)^2 - 2.1*x(1,:)^4 + 1/3*x(1,:)^6 + x(1,:)*x(2,:) - 4*x(2,:)^2 + 4*x(2,:)^4);
xMinA = [1; 2];
%Simulation parameters
eps = [10E-6, 10E-6]; maxIt = [1000, 1000];
%Grid and initial value
gridSize = 11; step = [1; 1]; 
x0 = cell(gridSize,gridSize); x0{1,1} = [-5; -5];
%All metrics
itRes = zeros(gridSize,gridSize,3); xRes = cell(gridSize,gridSize,3); 
fxRes = zeros(gridSize,gridSize,3);
trail = cell(gridSize,gridSize,3); 
 
for i=1:gridSize
    for j=1:gridSize
        %Filling grid of initial values
        x0{i,j} = x0{1,1} + [step(1)*(i-1); step(2)*(j-1)];
        %Launching algorithms and recording minima, iterations count and trail
        [xRes{i,j,1}, itRes(i,j,1), trail{i,j,1}] = Cauchy(f, x0{i,j}, eps, maxIt);
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
residual = -1*ones(gridSize, gridSize, 3);
minimaList = [];

for i=1:gridSize
    for j=1:gridSize
        for q=1:3
            hx = h(xRes{i,j,q});
            if round(norm(jF(xRes{i,j,q})), rC) == 0 && det(hx) > 0 && hx(1,1) > 0
                if xRes{i,j,q}(1) >= x0{1, 1}(1) && xRes{i,j,q}(2) >= x0{1, 1}(2) && xRes{i,j,q}(1) <= x0{end, end}(1) && xRes{i,j,q}(2) <= x0{end, end}(2)
                    exist = 0;
                    itRes(i,j,q) = abs(itRes(i,j,q));
                    for z=1:length(minimaList)
                        if round(norm(xRes{i,j,q} - minimaList(z).pos), ceil(rC/2)) == 0
                            residual(i,j,q) = norm(xRes{i,j,q} - minimaList(z).pos);
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
                        residual(i,j,q) = norm(xRes{i,j,q} - minimaList(end).pos);
                    end
                end
            else
                itRes(i,j,q) = NaN;
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
rightRes = [0; 0; 0]; wrongRes = [0; 0; 0];
convSpeed = ones(gridSize,gridSize,3);
genConvSpeed = [1; 1; 1]; genIt = [1; 1; 1];
genResidual = [0; 0; 0]; genElCount = [0; 0; 0]; 
wrongStepCount = [0; 0; 0];

for i=1:gridSize
    for j=1:gridSize
        for q=1:3
            if isnan(itRes(i,j,q))
                wrongRes(q) = wrongRes(q) + 1;
                convSpeed(i,j,q) = NaN;
            else
                elCount = 0;
                for p=3:length(trail{i,j,q})-1
                    if f(trail{i,j,q}{p}) - f(trail{i,j,q}{p-1}) < 0
                        coef = norm(f(trail{i,j,q}{p}) - f(trail{i,j,q}{p-1}))/ ...
                            norm(f(trail{i,j,q}{p+1}) - f(trail{i,j,q}{p}));
                        if coef ~= Inf && ~isnan(coef) && coef ~= 0 && coef < 10
                            convSpeed(i,j,q) = convSpeed(i,j,q) * coef;
                            elCount = elCount + 1;
                        end
                    else
                        if itRes(i,j,q) > 0
                            wrongStepCount(q) = wrongStepCount(q) + 1;
                        end
                    end
                end
                convSpeed(i,j,q) = convSpeed(i,j,q)^(1/(elCount));

                rightRes(q) = rightRes(q) + 1;
                if ~isempty(p)
                    genResidual(q) = genResidual(q) + residual(i,j,q);
                    genConvSpeed(q) = genConvSpeed(q) * convSpeed(i,j,q);
                    genIt(q) = genIt(q) * itRes(i,j,q);
                    genElCount(q) = genElCount(q) + 1;
                end
            end
        end
    end
end

for q=1:3
    genResidual(q) = genResidual(q) / genElCount(q);
    genConvSpeed(q) = genConvSpeed(q) ^ (1/genElCount(q));
    genIt(q) = genIt(q) ^ (1/genElCount(q));
end
%% Visualization for every starting point
%Plot parameters
plotGridFreq = 101;
thickeningCoef = 1.05;
lineThickness = 2;
for i=6
    for j=6
        X = x0{1,1}(1):0.1:x0{end,1}(1);
        Y = x0{1,1}(2):0.1:x0{1,end}(2);
        Z = fGrid(f, X, Y);
        L = nonLinspaceEnd(max(max(Z)), minimaList(minAIndex).f, 1.2, plotGridFreq);
        %Part of contour plot funcionality
        figure('Name', strcat("Contour", num2str(x0{i, j}(1)), num2str(x0{i, j}(2))));
        contourf(X, Y, Z, L, 'LineWidth', 0.5, 'ZLocation', minimaList(minAIndex).f - 1);
        clim([minimaList(minAIndex).f, 5]); hold on;
        plotTrail(trail{i,j,1},1, lineThickness);
        plotTrail(trail{i,j,2},2, lineThickness);
        plotTrail(trail{i,j,3},3, lineThickness);
        legend("","Cauchy","Newton","Markquardt");
        %Part of 3d plot funcionaluty
        figure('Name','3d');
        mesh(X, Y, Z); hold on;
        plotTrail3(f, trail{i,j,1},1, lineThickness);
        plotTrail3(f, trail{i,j,2},2, lineThickness);
        plotTrail3(f, trail{i,j,3},3, lineThickness);
        legend("","Cauchy","Newton","Markquardt");
    end
end
hold off;
%% Visualization of minimas
plotGridFreq = 51;
figure('Name',"All minimas");
hold on;
X = x0{1,1}(1):(x0{end,1}(1)-x0{1,1}(1))/100:x0{end,1}(1);
Y = x0{1,1}(2):(x0{1,end}(2)-x0{1,1}(2))/100:x0{1,end}(2);
Z = fGrid(f, X, Y);
L = nonLinspaceEnd(max(max(Z)), minimaList(minAIndex).f, 1.2, plotGridFreq);
for i = 1:length(minimaList)
    plot3(minimaList(i).pos(1), minimaList(i).pos(2), minimaList(i).f, '.','MarkerSize', 20, 'Color', [0.8, 0, 0]);
end
contourf(X, Y, Z, L, 'LineWidth', 0.5, 'ZLocation', minimaList(minAIndex).f - 1);
clim([minimaList(minAIndex).f, 5]);

for i = 1:length(minimaList)
    window = figure;   
    hold on;
    X = x0{1,1}(1):(x0{end,1}(1)-x0{1,1}(1))/100:x0{end,1}(1);
    Y = x0{1,1}(2):(x0{1,end}(2)-x0{1,1}(2))/100:x0{1,end}(2);
    Z = fGrid(f, X, Y);
    L = nonLinspaceEnd(max(max(Z)), minimaList(minAIndex).f, 1.2, plotGridFreq);
    contourf(X, Y, Z, L, 'LineWidth', 0.3, 'ZLocation', minimaList(minAIndex).f - 1);
    clim([minimaList(minAIndex).f, 5]);
    for j = 1:length(minimaList(i).trail)
        for z = 1:3
            if minimaList(i).trail(j).alg(z) == 1
                curTrail = trail{minimaList(i).trail(j).pos(1), minimaList(i).trail(j).pos(2), z};
                plotTrail(curTrail, z, 1.5);
            end
        end
    end
    plot(minimaList(i).pos(1), minimaList(i).pos(2), '.', 'MarkerSize', 20, 'Color', [1, 1, 1]);
    if i == minAIndex
        window.Name = 'Absolute';
    end
    hold off;
end
%% Colormap of expended function calculations
names = ["Cauchy", "Newton", "Markquard"];
X = zeros(gridSize,gridSize); Y = zeros(gridSize,gridSize);
for i=1:gridSize
    for j=1:gridSize
        X(i,j) = x0{i,j}(1);
        Y(i,j) = x0{i,j}(2);
    end
end
for i=1:3
    figure('Name',strcat(names(i)," Function calculation count"));
    pcolor(X,Y,itRes(:,:,i)); xlabel("x"); ylabel("y");
    colorbar
    figure('Name',strcat(names(i)," Average convergence speed"));
    pcolor(X,Y,convSpeed(:,:,i)); xlabel("x"); ylabel("y");
    colorbar
end











