%% Run grid test
f = @(x) (4*x(1,:)^2 - 2.1*x(1,:)^4 + 1/3*x(1,:)^6 + x(1,:)*x(2,:) - 4*x(2,:)^2 + 4*x(2,:)^4);
%Simulation parameters
eps = [10E-6, 10E-6]; maxIt = [1000, 1000];
step = 1; gridSize = 20; gridFreq = 51;
%All metrics
it = zeros(gridSize,gridSize,3); xMin = cell(gridSize,gridSize,3); fxMin = zeros(gridSize,gridSize,3);
trail = cell(gridSize,gridSize,3); 
%Grid and initial value
x0 = cell(gridSize,gridSize); x0{1,1} = [-10; -10];

for i=1:gridSize
    for j=1:gridSize
        %Filling grid of initial values
        x0{i,j} = x0{1,1} + [step*(i-1);step*(j-1)];
        %Launching algorithms and recording minima, iterations count and trail
        [xMin{i,j,1}, it(i,j,1), trail{i,j,1}] = Koshi(f, x0{i,j}, eps, maxIt);
        [xMin{i,j,2}, it(i,j,2), trail{i,j,2}] = Newton(f, x0{i,j}, eps(1), maxIt(1));
        [xMin{i,j,3}, it(i,j,3), trail{i,j,3}] = Mark(f, x0{i,j}, eps(1), maxIt(1));
        %Finding function value
        for q=1:3
            fxMin(i,j,q) = f(xMin{i,j,q});
        end
        %Visualization stuff
        if i == 3 && j == 3
            deltaX = x0{i,j}(1) - xMin{i,j,3}(1);
            deltaY = x0{i,j}(2) - xMin{i,j,3}(2);
            X = nonLinspaceMid(xMin{i,j,3}(1) - deltaX, xMin{i,j,3}(1) + deltaX, 1.1, gridFreq);
            Y = nonLinspaceMid(xMin{i,j,3}(2) - deltaY, xMin{i,j,3}(2) + deltaY, 1.1, gridFreq);
            Z = fGrid(f, X, Y);
            L = nonLinspaceEnd(max(max(Z)), fxMin(i,j,3), 1.2, 2*gridFreq);
            figure('Name','Contour'); 
            contour(X, Y, Z, L, 'LineWidth', 0.5); hold on;
            plotTrail(trail{i,j,1},'b-', 1.2);
            plotTrail(trail{i,j,2},'r-', 1.2);
            plotTrail(trail{i,j,3},'g-', 1.2);
            legend("","Koshi","Newton","Markquardt");

            figure('Name','3d'); 
            mesh(X, Y, Z); hold on;
            plotTrail3(f, trail{i,j,1},'b-', 1.2);
            plotTrail3(f, trail{i,j,2},'r-', 1.2);
            plotTrail3(f, trail{i,j,3},'g-', 1.2);
            legend("","Koshi","Newton","Markquardt");
        
        end
    end
end

hold off;
%% Results
rightResCount = 0;
X = zeros(gridSize,gridSize); Y = zeros(gridSize,gridSize);
for i=1:gridSize
    for j=1:gridSize
        X(i,j) = x0{i,j}(1);
        Y(i,j) = x0{i,j}(2);
    end
end
for i=1:3
    figure
    s = pcolor(X,Y,it(:,:,i)); xlabel("x"); ylabel("y");
    colorbar
end
for i=1:gridSize
    for j=1:gridSize
        for q=1:3
            if round(xMin{i,j,q}(1),3) == round(0.089842,3) && round(xMin{i,j,q}(2),3) == round(-0.712656,3)
                rightResCount = rightResCount + 1;
            end
        end
    end
end






