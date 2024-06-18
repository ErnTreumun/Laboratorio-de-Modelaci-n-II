function [x,val] = IntPointLinSolve(G,c,Adesg,bdesg,X0,max_iter)
% Minimizar 1/2*x^T*G*x + x^T*c
% sujeto a Adesg*x>= bdesg

    [m,n] = size(Adesg);
    x = X0;
    y = ones(m,1);
    l = ones(m,1);
    s = 0.1;

    for iter = 1:max_iter
        mu = 1/m * (y' * l);
        rd = G*x-(Adesg')*l+c;
        rp = Adesg * x - y - bdesg;
        rc = -diag(l) * y + s * mu * ones(m,1);

        R = [G, zeros(n,m), -(Adesg');
            Adesg, -eye(m), zeros(m,m);
            zeros(m,n), diag(l), diag(y)];

        r = [-rd; -rp; rc];
        sol = linsolve(R,r);

        
        dx = sol(1:n);
        dy = sol(n+1:n+m);
        dl = sol(n+m+1:end);

        % Busqueda lineal
        a = 1;
        while any(y+a*dy<=0) || any(l+a*dl<=0)
            a = a*0.9;
        end
        x = x + a*dx;
        y = y + a*dy;
        l = l + a*dl;
        
        % Valor objetivo
        val = (1/2) * x' * G * x + x' * c;
        % Imprimir el valor de la función objetivo en cada iteración (opcional)
        % fprintf('Iteración %d: Valor objetivo = %.6f\n', iter, val)

        s = 0.9*s;
    end
end