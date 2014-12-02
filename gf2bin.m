function X = gf2bin(x)
    X = zeros(x.m, size(x,2));
    x = x.x;
    for i = 1:size(X,1)
        X(i,:) = mod(x,2);
        x = (x - mod(x,2))/2;
    end
    X = gf(X);
end
