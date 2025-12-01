%% SIMULACIÓN DE LA ECUACIÓN DE VAN DER WAALS
% Este script consolida todas las simulaciones de la ecuación de van der Waals
% Autor: Angel Gabriel Pérez Sandoval
% Fecha: 2025

function van_der_waals_completo()
    clc;
    fprintf('=== SIMULACIÓN DE LA ECUACIÓN DE VAN DER WAALS ===\n\n');
    fprintf('Seleccione la simulación que desea ejecutar:\n');
    fprintf('1. Isotermas básicas (t = 0.8 a 1.2)\n');
    fprintf('2. Isoterma con regla de Maxwell (requiere temperatura)\n');
    fprintf('3. Diagrama completo de fases con curva de coexistencia\n');
    fprintf('4. Ejecutar todas las simulaciones\n');
    fprintf('0. Salir\n\n');
    
    opcion = input('Ingrese su opción: ');
    
    switch opcion
        case 1
            simulacion_isotermas_basicas();
        case 2
            simulacion_maxwell();
        case 3
            simulacion_diagrama_fases();
        case 4
            simulacion_isotermas_basicas();
            fprintf('\nPresione Enter para continuar...\n');
            pause;
            simulacion_maxwell();
            fprintf('\nPresione Enter para continuar...\n');
            pause;
            simulacion_diagrama_fases();
        case 0
            fprintf('Simulación finalizada.\n');
            return;
        otherwise
            fprintf('Opción no válida.\n');
    end
end

%% SIMULACIÓN 1: Isotermas básicas
function simulacion_isotermas_basicas()
    fprintf('\n--- SIMULACIÓN 1: Isotermas básicas ---\n');
    
    figure('Name', 'Isotermas de van der Waals', 'NumberTitle', 'off');
    f = @(t,x) 8*t./(3*x-1)-3./(x.^2);
    hold on
    axis([0.5 4 0 2])
    
    colores = ['b', 'g', 'k', 'm', 'c'];
    idx = 1;
    
    for t = 0.8:0.1:1.2
        v = linspace(0.5, 4, 100);
        p = f(t, v);
        plot(v, p, 'Color', colores(idx), 'LineWidth', 1.5, ...
             'DisplayName', sprintf('t = %.1f', t));
        idx = idx + 1;
    end
    
    plot(1, 1, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r', ...
         'DisplayName', 'Punto crítico');
    xlabel('V_R (Volumen reducido)', 'FontSize', 12)
    ylabel('P_R (Presión reducida)', 'FontSize', 12)
    title('Isotermas de la ecuación de van der Waals', 'FontSize', 14)
    legend('Location', 'northeast')
    grid on
    hold off
    
    fprintf('Gráfica generada: Isotermas para t = 0.8, 0.9, 1.0, 1.1, 1.2\n');
    fprintf('El punto rojo marca el punto crítico (V_R = 1, P_R = 1, T_R = 1)\n');
end

%% SIMULACIÓN 2: Isoterma con regla de Maxwell
function simulacion_maxwell()
    fprintf('\n--- SIMULACIÓN 2: Isoterma con regla de Maxwell ---\n');
    
    t = input('Temperatura reducida (0.8-0.99): ');
    
    if t >= 1 || t < 0.8
        fprintf('Temperatura fuera de rango. Use valores entre 0.8 y 0.99\n');
        return;
    end
    
    figure('Name', 'Regla de Maxwell', 'NumberTitle', 'off');
    
    % Gráfica de la isoterma
    f = @(x) 8*t./(3*x-1)-3./(x.^2);
    Vr = linspace(0.5, 4, 100);
    Pr = f(Vr);
    hold on
    axis([0 4 0 1.4])
    plot(Vr, Pr, 'r', 'LineWidth', 2, 'DisplayName', sprintf('Isoterma t = %.2f', t))
    xlabel('V_R (Volumen reducido)', 'FontSize', 12)
    ylabel('P_R (Presión reducida)', 'FontSize', 12)
    
    % Cálculo de máximo y mínimo
    pol = [1 -9/(4*t) 3/(2*t) -1/(4*t)];
    z = sort(raices_3(pol));
    plot(z(1), f(z(1)), 'k*', 'MarkerSize', 10, 'DisplayName', 'Extremos locales')
    plot(z(2), f(z(2)), 'k*', z(3), f(z(3)), 'k*', 'MarkerSize', 10, 'HandleVisibility', 'off')
    
    if f(z(2)) > 0
        a = f(z(2)) + 0.01;
    else
        a = 0.01;
    end
    b = f(z(3)) - 0.01;
    
    % Calcula p para que las áreas sean iguales
    f1 = @(p) igualArea(p, t);
    p = fzero(f1, [a b]);
    
    plot([0 4], [p p], 'k--', 'LineWidth', 1.5, 'DisplayName', 'Presión de equilibrio')
    pol = [1 -(p+8*t)/(3*p) 3/p -1/p];
    v = sort(raices_3(pol));
    plot(v(1), p, 'bo', v(2), p, 'bo', v(3), p, 'bo', ...
         'MarkerSize', 8, 'MarkerFaceColor', 'b', 'DisplayName', 'Puntos de corte')
    
    % Rellenar las áreas
    xx = [v(1) Vr(Vr >= v(1) & Vr <= v(2)) v(2)];
    yy = [f(v(1)) Pr(Vr >= v(1) & Vr <= v(2)) f(v(2))];
    fill(xx, yy, 'y', 'FaceAlpha', 0.5, 'EdgeColor', 'none', ...
         'DisplayName', 'Áreas iguales (Maxwell)')
    
    xx = [v(2) Vr(Vr >= v(2) & Vr <= v(3)) v(3)];
    yy = [f(v(2)) Pr(Vr >= v(2) & Vr <= v(3)) f(v(3))];
    fill(xx, yy, 'y', 'FaceAlpha', 0.5, 'EdgeColor', 'none', ...
         'HandleVisibility', 'off')
    
    plot(Vr, Pr, 'r', 'LineWidth', 2, 'HandleVisibility', 'off')
    
    % Imprimir datos
    text(1, 1.35, sprintf('P_R = %.4f', p), 'FontSize', 10, 'FontWeight', 'bold')
    text(1, 1.25, sprintf('V_{gas} = %.4f', v(3)), 'FontSize', 10)
    text(1, 1.15, sprintf('V_{líquido} = %.4f', v(1)), 'FontSize', 10)
    
    title('Regla de Maxwell - Equilibrio líquido-gas', 'FontSize', 14)
    legend('Location', 'northeast')
    grid on
    hold off
    
    fprintf('\nResultados:\n');
    fprintf('  Presión de equilibrio: %.4f\n', p);
    fprintf('  Volumen líquido: %.4f\n', v(1));
    fprintf('  Volumen gas: %.4f\n', v(3));
    fprintf('  Las áreas amarillas son iguales (Regla de Maxwell)\n');
end

%% SIMULACIÓN 3: Diagrama completo de fases
function simulacion_diagrama_fases()
    fprintf('\n--- SIMULACIÓN 3: Diagrama de fases completo ---\n');
    
    figure('Name', 'Diagrama de Fases', 'NumberTitle', 'off');
    Vr = linspace(0.1, 4, 100);
    hold on
    i = 1;
    PV = zeros(9, 2);
    axis([0.4 4 0.2 1.4])
    
    for t = 0.8:0.05:1.2
        f = @(x) 8*t./(3*x-1)-3./(x.^2);
        
        if t < 1
            % Calcula máximo y mínimo local de p
            pol = [1 -9/(4*t) 3/(2*t) -1/(4*t)];
            z = sort(raices_3(pol));
            
            if f(z(2)) > 0
                a = f(z(2)) + 0.01;
            else
                a = 0.01;
            end
            b = f(z(3)) - 0.01;
            
            % Calcula p para que las áreas sean iguales
            f1 = @(p) igualArea(p, t);
            p = fzero(f1, [a b]);
            pol = [1 -(p+8*t)/(3*p) 3/p -1/p];
            v = sort(raices_3(pol));
            
            % Guarda Vlíquido, Vgas y la presión constante p
            PV(i, 1) = v(1); PV(i, 2) = p;
            PV(i+1, 1) = v(3); PV(i+1, 2) = p;
            i = i + 2;
            
            % Gráfica
            Pr = f(Vr);
            Vf = [Vr(Vr <= v(1)) v(1) v(3) Vr(Vr >= v(3))];
            Pf = [Pr(Vr <= v(1)) f(v(1)) f(v(3)) Pr(Vr >= v(3))];
            plot(Vf, Pf, 'r', 'LineWidth', 1.2)
        else % En el punto crítico y por encima
            Pr = f(Vr);
            plot(Vr, Pr, 'r', 'LineWidth', 1.2)
        end
    end
    
    % Une los puntos (p Vlíquido), (p Vgas) - cambio de fase
    PV(i, 1) = 1; PV(i, 2) = 1; % Punto crítico
    PVs = sortrows(PV, 1);
    plot(PVs(:, 1), PVs(:, 2), 'b', 'LineWidth', 2.5, 'DisplayName', 'Curva de coexistencia')
    plot(1, 1, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b', ...
         'DisplayName', 'Punto crítico')
    
    % Anotar regiones
    text(0.7, 0.8, 'LÍQUIDO', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'blue')
    text(2.5, 0.8, 'GAS', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'blue')
    text(1.5, 0.5, 'COEXISTENCIA', 'FontSize', 10, 'Color', 'blue', 'Rotation', 30)
    
    xlabel('V_R (Volumen reducido)', 'FontSize', 12)
    ylabel('P_R (Presión reducida)', 'FontSize', 12)
    title('Diagrama P-V con curva de coexistencia líquido-gas', 'FontSize', 14)
    legend('Location', 'northeast')
    grid on
    hold off
    
    fprintf('Gráfica generada: Diagrama completo de fases\n');
    fprintf('  - Isotermas rojas: t = 0.8 a 1.2\n');
    fprintf('  - Curva azul: frontera de coexistencia líquido-gas\n');
    fprintf('  - Punto crítico: (V_R = 1, P_R = 1, T_R = 1)\n');
end

%% FUNCIÓN: Raíces de ecuación cúbica
function x = raices_3(p)
    Q = (p(2)*p(2) - 3*p(3))/9;
    R = (2*p(2)^3 - 9*p(2)*p(3) + 27*p(4))/54;
    x = zeros(3, 1);
    
    if (R*R) < (Q^3)
        tetha = acos(R/sqrt(Q^3));
        x(1) = -2*sqrt(Q)*cos(tetha/3) - p(2)/3;
        x(2) = -2*sqrt(Q)*cos((tetha+2*pi)/3) - p(2)/3;
        x(3) = -2*sqrt(Q)*cos((tetha-2*pi)/3) - p(2)/3;
    else
        A = -sign(R)*nthroot(abs(R) + sqrt(R*R - Q^3), 3);
        if A == 0
            B = 0;
        else
            B = Q/A;
        end
        x(1) = (A+B) - p(2)/3;
        x(2) = -(A+B)/2 - p(2)/3 + (sqrt(3)*(A-B)/2)*sqrt(-1);
        x(3) = -(A+B)/2 - p(2)/3 - (sqrt(3)*(A-B)/2)*sqrt(-1);
    end
end

%% FUNCIÓN: Igualar áreas (Regla de Maxwell)
function res = igualArea(p, t)
    pol = [1 -(p+8*t)/(3*p) 3/p -1/p];
    v = sort(raices_3(pol));
    
    area_1 = p*(v(2)-v(1)) - 8/3*t*log((3*v(2)-1)/(3*v(1)-1)) - ...
             3*(v(1)-v(2))/(v(1)*v(2));
    area_2 = 8/3*t*log((3*v(3)-1)/(3*v(2)-1)) + 3*(v(2)-v(3))/(v(2)*v(3)) - ...
             p*(v(3)-v(2));
    
    res = area_1 - area_2;
end
