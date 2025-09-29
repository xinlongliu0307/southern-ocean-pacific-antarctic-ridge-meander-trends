function slope = Sen_Slope(x)
x = x(:); n = length(x);
slopes = [];
for i = 1:n-1
    for j = i+1:n
        slopes(end+1) = (x(j) - x(i)) / (j - i);
    end
end
slope = median(slopes, 'omitnan');
end
