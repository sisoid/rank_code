function Y = coder(X)

%определяем параметры кода
[n, k] = size(X);
d = n - k + 1;

%преобразуем информационную матрицу в вектор
x = gf(2.^(0:n-1), n)*X.x;

%задаем проверочную матрицу
H = gf(zeros(d-1,n), n);
%первая строка проверочной матрицы
H(1,:) = gf(2.^(0:n-1), n);
%остальные строки - ее фробениусовские степени
for i = 2:d-1
    H(i,:) = H(i-1,:).^2;
end

%генерируем порождаующую матрицу в систематическом виде
G = gf(eye(k), n);
%должно выполняться условие G*H'=0
for i = 1:k
    G(i,k+1:n) = (H(:,k+1:n)^-1)*H(:,i);
end

%генерируем кодовый вектор
y = x*G;

%преобразуем кодовый вектор в матричную форму
Y = gf2bin(y);

end
