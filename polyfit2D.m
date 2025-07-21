function result = polyfit2D(phasec, mask, n)
    [y, x] = find(mask);
    
    A = zeros(length(x), (n+1)*(n+2)/2);
    col = 0;
    for i = 1:n
        for j = 1:i
            col = col + 1;
            A(:, col) = x.^(i-j) .* y.^j;
        end
    end
    
    result = inv(A'*A)*A'*phasec(mask==1);
end
