function res = nonLinspaceMid(min, max, c, n)
    if mod(n,2) == 1
        x = nonLinspace((max - min)/2, c, n/2);
        xl = min + x;
        xr = max - wrev(x);
        res = [xl xr(2:end)];
    else
        x = nonLinspace((max - min)/2, c, n/2 - 1);
        xl = min + x;
        xr = max - wrev(x);
        res = [xl xr(2:end)];
    end
end

