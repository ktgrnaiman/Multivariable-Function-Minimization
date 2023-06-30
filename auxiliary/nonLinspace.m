function res = nonLinspace(d, c, n)
    x = zeros(1,n); sum = 0;
    for i=0:n-2
        sum = sum + c^(-i);
    end
    x(2) = d/sum;
    for i=3:n
        x(i) = x(i-1) + (x(i-1) - x(i-2))/c;
    end
    res = x;
end
