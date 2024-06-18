function [x, objval] = FWOptStep(Q, c, Adesg, bdesg, Aeq, beq, xmin, xmax, max_iter, X0,Opt)
    % Minimizando una función de la forma
    % x^T Q x + c^T x sujeto a
    % Adesg x <= bdesg; Aeq x = beq; xmín <= x <= xmáx
    % xmin, Aquí Opt permite utilizar o no un paso óptimo
    %Con opt= 0 NO hacemos paso óptimo y con Opt = 1 sí
    % <= x <=xmax partiendo de X0 y con un máximo de iteracioens max_iter
    % Inicialización
    [~, n] = size(Q);
    x = X0;

    % Opciones para suprimir mensajes de salida de linprog
    options = optimoptions('linprog', 'Display', 'off');
    
    for iter = 1:max_iter
        % Calcular el gradiente
        gradient = (Q' +Q) * x + c;
        
        % Resolver subproblema de minimización restringida (búsqueda lineal)
        direction= linprog(gradient,Adesg,bdesg,Aeq,beq,xmin,xmax,options);
        if Opt==0
            % Calcular el tamaño del paso (paso óptimo)
            step_size = 2 / (iter + 2);

        else
            A1 = x' * Q * x;
            B1 = x'*(Q + Q')*direction;
            C1 = direction'*Q*direction;
            alpha = A1 - B1 + C1;
            beta = B1 - 2*A1 + c'*(direction - x);

            step_size = QP1D(alpha,beta,0, 0, 1);
        end

        % Actualizar la solución
        x = (1- step_size)*x + step_size * direction;
        
        % Calcular el valor objetivo y la norma del gradiente
        objval = x'*Q*x + c'*x;
    end
end
