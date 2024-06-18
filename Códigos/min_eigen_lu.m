function lambda_min = min_eigen_lu(P, L, U, num_iter, tol)
    
    n = size(P, 1);
    x = randn(n, 1);  % Vector inicial aleatorio
    x = x / norm(x);  % Normalizar el vector inicial
    
    for k = 1:num_iter

		z = lu_solve(P,L,U,x);

        % Normalizar y encontrar el nuevo vector de iteración
        mu = norm(z);
        x_next = z / mu;
        
        % Verificar la convergencia
        if norm(x_next - x) < tol
            break;
        end
        
        x = x_next;
    end
    
    % El valor propio mínimo es el inverso de mu
    lambda_min = 1 / mu;
end