function Z = fGrid(f, x, y)
    Z = zeros(size(x,2), size(y,2));
    for i=1:size(x,2)
        for j=1:size(y,2)
            Z(i, j) = f([x(i); y(j)]);
        end
    end
    Z;
end