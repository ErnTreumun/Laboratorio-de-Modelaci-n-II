function [x, objval] = QP1d(a, b, c, x0, x1)
    % Esta función resuelve un problema de minimización cuadrático
    % unidimensional con límites x0 y x1, donde la función objetivo es:
    % f(x) = a*x^2 + b*x + c

    if a <= 0
        if a*x0^2 + b*x0 < a*x1^2 + b*x1
            x = x0;
        else
            x = x1;
        end
    else
        x_vertex = -b / (2 * a);

        if x0 < x_vertex && x_vertex < x1
            x = x_vertex;
        else
            if a*x0^2 + b*x0 < a*x1^2 + b*x1
                x = x0;
            else
                x = x1;
            end
        end
    end

    objval = a * x^2 + b * x + c;
end
