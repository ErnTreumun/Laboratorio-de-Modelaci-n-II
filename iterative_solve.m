function [x,cont] = iterative_solve(G, P,L,U, A, l, y, b, tol, max_iter, lambdamingG)
	% Resolver (G + A^\top diag(y)^{-1} diag(l) A)x = b + (términos que dependen de y).
    % Aquí PG = LU.
	% b: Vector de términos independientes
    % x0: Vector de estimación inicial, tol: Tolerancia para la convergencia
    % max_iter: Número máximo de iteraciones
	% A,l,y,b son los asociados a nuestro problema de interior point
	
    [m,n] = size(A);

    % Extraer bloques de b
    b1 = b(1:n);
    b2 = b(n+1:n+m);
    b3 = b(n+m+1:end);

	%Esto es necesario para verificar si el método de Jacobi converge
	D = (1./y).*l;
	norm1A2= norm(A,1)^2;
	%Si el radio espectral es menor que 1, el método iterativo sí converge.

	if max(abs(D))*m*norm1A2*lambdamingG < 1

		x1 = randn(n, 1);  % Vector inicial aleatorio
		iter = 0;

	    while iter < max_iter
        	x1_old = x1;
        	%Resolver para x2 dado x1
        	x2 = A*x1_old - b2;
	
        	% Resolver para x3 dado x2
        	x3 = (b3 - l.*x2) ./ y;
        	
        	% Resolver para x1 dado x3
        	y1 = (A')*x3 + b1;
        	
        	x1 = lu_solve(P,L,U,y1);
	
        	% Verificar la convergencia
        	if norm(x1 - x1_old, inf) < tol
				x= [x1; x2; x3];
            	return;
        	end
        
        iter = iter + 1;

		end

		return
	end
	%max(abs(D))*norm1A*lambdamingG
	%Quiero resolver Gx1 = - Aaux x1 + baux

	%baux = A^\top diag(y)^{-1} (b3 + diag(l) b2) + b1;
	baux= A'*( (b3 + l.*b2)./y) + b1;
	%Aaux = A^\top (diag(y^{-1}) diag(l) A)
	Aaux = A'* diag((1 ./ y) .* l )*A;

	[P2,L2,U2] = lu_decomposition(G+ Aaux);
	x1 = lu_solve(P2,L2,U2, baux);

	x2 = A*x1 - b2;
	x3 = (b3 - l.*x2)./y;
	x= [x1; x2; x3];
end