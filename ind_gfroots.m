function roots = ind_gfroots(coeffs, free)
    %размерность поля
    n = coeffs.m;
    
    %размерность многочлена
    t = size(coeffs, 2);
    
    %вектор переменных
    xvar = gf(2,n).^(1:2^n-1);
    
    %считаем многочлен для каждой переменной
    expr = xvar.^(2^t) + free;
    for i=1:t
        expr = expr + coeffs(i)*xvar.^(2^(i-1));
    end
    
    %если уровнение имеет свободный член, возвращаем первый корень
    if free ~= 0
        roots = gf(2, n)^find(expr == 0, 1);
        return;
    end
    
    %находим все корни
    all_roots = gf2bin(gf(2, n).^find(expr == 0));
    
    %находим линейно независимые корни    
    i = 1;
    j = 1;
    all_r = rank(all_roots);
    ind_roots = gf(zeros(n, all_r));
    r = 0;
    while r < all_r
        ind_roots(:,i) = all_roots(:,j);
        if rank(ind_roots) > r
            i = i + 1;
            r = r + 1;
        end
        j = j + 1;
    end
    
    roots = gf(2.^(0:n-1), n)*ind_roots.x;
end
