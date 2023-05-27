function plotTrail(trail, LineProperties, LineWidth)
x = zeros(1,length(trail)); y = x;
for i=1:length(trail)
    x(i) = trail{i}(1); y(i) = trail{i}(2);
end
plot(x,y,LineProperties, 'LineWidth', LineWidth);
end