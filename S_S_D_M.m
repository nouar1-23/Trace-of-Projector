function [D1, l1, rrr] =S_S_D_M(a, b)
% S_S_D_M Анализ спектральных свойств матрицы на основе разложения Шура.
%
% Входные параметры:
%   a, b - Диапазон изменения параметра сканирования
%
% Выходные параметры:
%   D1 - Массивы значений параметра 
%   l1 - Массивы невязок (след проектора)
%   rrr - Значения параметров, при которых изменяется размерность подпространства

    [A_original] = matre_x(); % Получение исходной матрицы

    % --- Глобальные параметры вычислений ---
    h = 0.01;           % Шаг сканирования по параметру λ
    w_max = 1e25;       % Пороговое ограничение нормы решения
    u_max = 1e20;       % Ограничение числа обусловленности
    ip = 1e-15;         % Машинная точность (критерий сходимости)
    
    [n, ~] = size(A_original);
    I = eye(n);
    A_qz = A_original;
    
    r1_val = a; 
    i = 1; 
    u = 1;
    last_k_dim = -1; 
    rrr = [];

    % Статическое разложение Шура (выполняется один раз для оптимизации)
    [Z_fixed, T_fixed] = schur(A_qz, 'complex'); 

    % --- Основной цикл сканирования ---
    while r1_val <= b
        % Нормировка текущей матрицы
        AA_current = T_fixed / r1_val;
        
        % Выбор собственных значений внутри единичного круга
        select = abs(diag(AA_current)) < 1.0;
        current_k_dim = sum(select);
        
        % Отслеживание точек изменения размерности инвариантного подпространства
        if current_k_dim ~= last_k_dim
            rrr(u) = r1_val;
            u = u + 1;
        % Пересчет проектора и матрицы Сильвестера при необходимости
            [Z_s, AA_s] = ordschur(Z_fixed, AA_current, select);
            k_dim = sum(select);
            q_dim = n - k_dim;
            last_k_dim = k_dim;
            
            A11 = AA_s(1:k_dim, 1:k_dim) * r1_val;
            A12 = AA_s(1:k_dim, k_dim+1:n) * r1_val;
            A22 = AA_s(k_dim+1:n, k_dim+1:n) * r1_val;
            
            % Решение матричного уравнения Сильвестера
            M_sys = sylvester(A11/r1_val, -A22/r1_val, -A12/r1_val);

            I_k = eye(k_dim);
            I_q = eye(q_dim);
            
            % Матрица преобразования
            inv_Q1_transform = [I_k, -M_sys; zeros(q_dim, k_dim), I_q] * Z_s';
            
            % Контроль обусловленности блока A22
            if cond(A22) > u_max
                m1(i) = NaN; l1(i) = NaN;
                e1(i) = r1_val;
                disp(['Предупреждение: Точка r = ', num2str(r1_val), ' критически близка к спектру']);
                r1_val = r1_val + h;
                i = i + 1;
                continue 
            end 
            
            w22_sub_qz = M_sys' * M_sys + I_q;
            inv_BB = I_q / A22; 
            w11_sub_qz = I_k;
            
            Z1 = Z_s(:, 1:k_dim); 
            Z2 = Z_s(:, k_dim+1:n);
            % Вычисление спектрального проектора P
            pr = Z1 * (Z1' - M_sys * Z2');
        end

        % --- Итерационный процесс для блока X (метод удвоения) ---
        aa_iter = A11 / r1_val;
        x0_mat = w11_sub_qz;
        x1_mat = x0_mat + (aa_iter') * x0_mat * aa_iter;
        s_count = 0;
        
        while (k_dim > 0 && max(abs(x1_mat - x0_mat), [], 'all') > ip && ...
               max(abs(x1_mat), [], 'all') < w_max && s_count < 100)
            s_count = s_count + 1;
            for j = 1:3
                aa_iter = aa_iter * aa_iter;
                x0_mat = x1_mat;
                x1_mat = x0_mat + aa_iter' * x0_mat * aa_iter;
            end
        end

        % --- Итерационный процесс для блока Y ---
        bb_iter = r1_val * inv_BB;
        y0_mat = bb_iter' * w22_sub_qz * bb_iter;
        y1_mat = y0_mat + bb_iter' * y0_mat * bb_iter;
        s_count = 0;
        
        while (k_dim < n && max(abs(y1_mat - y0_mat), [], 'all') > ip && ...
               max(abs(y1_mat), [], 'all') < w_max && s_count < 100)
            s_count = s_count + 1;
            for j = 1:3
                bb_iter = bb_iter * bb_iter;
                y0_mat = y1_mat;
                y1_mat = y0_mat + bb_iter' * y0_mat * bb_iter;
            end
        end
        
        % --- Формирование итоговой матрицы Грамма H ---
        mm_zero = zeros(k_dim, n - k_dim);
        zz_block = [x1_mat, mm_zero; mm_zero', y1_mat];
        H1_final = inv_Q1_transform' * zz_block * inv_Q1_transform;
        
        % Сбор статистических данных
        m1(i) = max(abs(H1_final), [], 'all');
        l1(i) = trace(pr); % След проектора (соответствует k_dim)

        e1(i) = r1_val;
        i = i + 1;
        r1_val = r1_val + h;
    end
    
    D1 = e1;
