function res = nonLinspaceMid(min, max, c, n)
    if mod(n,2) == 1
        x = nonLinspace((max - min)/2, c, (n+1)/2);
        xl = min + x;
        xr = max - wrev(x);
        res = [xl xr(2:end)];
    else
        error('Vector size must be odd number');
    end
end

