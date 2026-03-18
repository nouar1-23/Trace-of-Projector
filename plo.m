function [] = plo()
    % --- Выполнение расчета ---
    % Параметры: сканирование от r = 0.9 до r = 2
    [F1, l1, rrr] = S_S_D_M(0.9, 2);

    rrr=rrr(2:end);

    % --- Построение графика размерности (Trace of Projector) ---
    figure 
    hold on;

    % Отрисовка основной линии следа проектора
    plot(F1, l1, 'b-', 'LineWidth', 4); 

    % Цикл для нанесения меток собственных значений (точек смены размерности)
    for i = 1:length(rrr)
        x_val = rrr(i);
        
        % Интерполяция для нахождения значения y в точке перехода
        try
            y_val = interp1(F1, l1, x_val, 'previous');
        catch
            y_val = 0; 
        end
        
        % Отрисовка точки (собственного значения) на оси абсцисс
        plot(x_val, 0, 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'Color', 'k');
        
        % Отрисовка вертикальной пунктирной линии связи
        line([x_val x_val], [0 y_val], 'Color', 'k', 'LineStyle', ':', 'LineWidth', 3);
        
        % Добавление текстовой метки значения r (с поворотом для читаемости)
        text(x_val + 0.01, 0.005 * max(l1), num2str(x_val, '%.2f'), ...
            'Color', 'k', 'FontSize', 20, 'Rotation', 45);
    end

    % --- Оформление осей и шрифтов ---
    xlabel('r', 'FontSize', 25);
    ylabel('tr(P)', 'Interpreter', 'latex');
    
    % Установка размера шрифта для делений осей
    set(gca, 'FontSize', 25);
        
end