function Y = decoder(Z, d, A, C)

%определяем параметры кода
n = size(Z, 1);

%задаем проверочную матрицу
H = gf(zeros(d-1,n), n);
%первая строка проверочной матрицы
H(1,:) = gf(2.^(0:n-1), n);
%остальные строки - ее фробениусовские степени
for i = 2:d-1
    H(i,:) = H(i-1,:).^2;
end

%преобразуем кодовую матрицу в вектор
z = gf(2.^(0:n-1), n)*Z.x;

%вычисляем синдром
s = z*H';

%количество стираний столбцов
l = size(C, 1);

%количество стираний строк
v = size(A, 2);

%если есть столбцевые стирания
if l > 0  
    %находим матрицу маленьких g
    g = gf(2.^(0:n-1), n)*C.x';
    gg = gf(zeros(l), n);
    for i = 1:l
        gg(i,:) = g.^(2^(i-1));
    end
    
    %находим гаммы
    Gamma = gg(1, :).^(2^l)/gg;
    Gamma(l+1) = gf(1, n);
    
    %находим матрицу из гамм
    G = gf(zeros(d-1, d-1-l), n);
    for i = 1:d-1
        for j = 1:d-1-l
            if 0 <= i-j && i-j <= l
                G(i,j) = Gamma(i-j+1)^(2^(j-1));
            else
                G(i,j) = gf(0, n);
            end
        end
    end
      
    %модифицируем синдром
    s = s*G;
end

if v > 0   
    %модифицируем синдром
    s = s.^(2.^(n:-1:n-(d-1-l)+1));
    s(1:end) = s(end:-1:1);
    
    %находим матрицу маленьких a
    a = (gf(2.^(0:n-1), n)*A.x).^(2.^(n-(d-1-l)+1));
    aa = gf(zeros(v), n);
    for i = 1:v
        aa(i,:) = a.^(2^(i-1));
    end
    
    %находим лямбды
    Lambda = aa(1, :).^(2^v)/aa;
    Lambda(v+1) = gf(1, n);

    %находим матрицу из лямбд
    L = gf(zeros(d-1-l, d-1-l-v), n);
    for i = 1:d-1-l
        for j = 1:d-1-l-v
            if 0 <= i-j && i-j <= v
                L(i,j) = Lambda(i-j+1)^(2^(j-1));
            else
                L(i,j) = gf(0, n);
            end
        end
    end
       
    %модифицируем синдром
    s = s*L;
end

%если синдром ненулевой, приступаем к декодированию рангового кода
if rank(s) > 0   
    %начальное значение параметра t
    t = floor((d-1-l-v)/2);
    
    %начальный вид матрицы M
    M = gf(zeros(t), n);
    for i=1:t
        for j=1:t
            M(i,j) = s(i+j-1)^(2^(n-j+1));
        end
    end
    
    %уменьшаем размер матрицы M, пока она не станет невырожденной
    while det(M) == 0 && t > 1
        t = t - 1;
        M = M(1:t, 1:t);
    end
    
    %если матрица M всегда вырождена - отказ от декодирования
    if det(M) == 0
        Y = Z;
        return;
    end
    
    %находим сигмы
    sigma = s(t+1:2*t).^(2.^(n:-1:n-t+1)) / M;
    
    %находим линейно независимые корни
    locators = ind_gfroots(sigma, 0);
    
    %если число корней меньше необходимого - отказ от декодирования
    if size(locators, 2) < t
        Y = Z;
        return;
    end
    
    %находим базис ошибок
    XX=gf(zeros(t), n);
    for i=1:t
        XX(:,i) = locators'.^(2^(i-1));
    end
    basis = s(1:t)/XX;
    
    %если были стирания строк
    if v > 0
        %меняем местами локаторы и базис ошибок
        locators = locators + basis;
        basis = locators - basis;
        locators = locators - basis;
        
        %убираем влияние матрицы из лямбд
        for i = 1:t
            all_roots = ind_gfroots(Lambda(1:v), basis(i));
            basis(i) = all_roots(1);
        end
        
        %убираем степень
        basis = basis.^(2^((d-1-l)-1));
    end
    
    %если были стирания столбцов
    if l > 0
        %убираем влияние матрицы из гамм
        for i = 1:t
            all_roots = ind_gfroots(Gamma(1:l), locators(i));
            locators(i) = all_roots(1);
        end
    end
    
    Z = Z - gf2bin(basis)*gf2bin(locators)';
end

%тут будет исправление ранговой ошибки

%исправляем стирания строк
if v > 0
    %преобразуем кодовую матрицу в вектор
    z = gf(2.^(0:n-1), n)*Z.x;
    
    %вычисляем синдром
    s = z*H';
    
    %исключаем стирания столбцов, если они были
    if l > 0
        s = s*G;
    end
    
    %модифицируем синдром
    s = s.^(2.^(n:-1:n-(d-1-l)+1));
    s(1:end) = s(end:-1:1);
    
    %исправляем стирания строк
    theta = s(1:v)/aa';
    if l > 0
        for i = 1:v
            all_roots = ind_gfroots(Gamma(1:l), theta(i));
            theta(i) = all_roots(1);
        end
    end
    R = gf2bin(theta)';
    Z = Z - A*R;
end

%исправляем стирания столбцов
if l > 0
    %преобразуем кодовую матрицу в вектор
    z = gf(2.^(0:n-1), n)*Z.x;
    
    %вычисляем синдром
    s = z*H';
        
    %исправляем стирания столбцов
    w = s(1:l)/gg';
    W = gf2bin(w);
    Z = Z - W*C;
end

Y = Z;
