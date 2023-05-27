function plotTrail3(f, trail, LineProperties, LineWidth)
x = zeros(1,length(trail)); y = x; z = x;
for i=1:length(trail)
    x(i) = trail{i}(1); y(i) = trail{i}(2);
    z(i) = f(trail{i});
end
plot3(x,y,z,LineProperties, 'LineWidth', LineWidth);
end