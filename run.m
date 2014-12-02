%параметры кода
n = 10;
k = 2;
d = n - k + 1;

Err = 0;
for it = 1:10
    %генерируем информационный вектор
    X = gf(randi(2, n, k) - 1);
    
    %кодируем
    Y = coder(X);
    
    %генерируем ранговые ошибки
    err_rank = randi(floor((d-1)/2)+1)-1;
    E = gf(zeros(n, n));
    while rank(E) ~= err_rank
        E = gf(randi(2, n, err_rank) - 1)*gf(randi(2, err_rank, n) - 1);
    end
    
    %генерируем стирания строк
    err_row = randi(d - 2*err_rank)-1;
    A = gf(zeros(n, err_row));
    R = gf(zeros(err_row, n));
    while rank(A) ~= err_row
        A = gf(randi(2, n, err_row) - 1);
        R = gf(randi(2, err_row, n) - 1);
    end
    
    %генерируем стирания столбцов
    err_column = randi(d - 2*err_rank - err_row)-1;
    W = gf(zeros(n, err_column));
    C = gf(zeros(err_column, n));
    while rank(C) ~= err_column
        W = gf(randi(2, n, err_column) - 1);
        C = gf(randi(2, err_column, n) - 1);
    end
    
    %формируем искаженное сообщение
    Z = Y + A*R + W*C + E;
    
    %[err_row err_column err_rank err_row+err_column+2*err_rank]
    Err = Err + rank(Y - decoder(Z, d, A, C));
end

Err
