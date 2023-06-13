function Z = fGrid(f, x, y)
    Z = zeros(length(x), length(y));
    for i=1:length(y)
        for j=1:length(x)
            Z(i, j) = f([x(j); y(i)]);
        end
    end
    Z;
end

