function x = lu_solve(P,L,U,b)
    % x = lu_solve(A, b) is the solution to A x = b
    % A must be a square matrix
    % b must be a vector of the same leading dimension as L
    z = P*b;
    y = forward_sub(L, z);
    x = back_sub(U, y);
end

% Función back_sub
function x = back_sub(U, y)
    % x = back_sub(U, y) is the solution to U x = y
    % U must be an upper-triangular matrix
    % y must be a vector of the same leading dimension as U

    n = size(U, 1);
    x = zeros(n, 1);
    for i = n:-1:1
        tmp = y(i);
        for j = i+1:n
            tmp = tmp - U(i,j) * x(j);
        end
        x(i) = tmp / U(i,i);
    end
end

% Función forward_sub
function y = forward_sub(L, b)
    % y = forward_sub(L, b) is the solution to L y = b
    % L must be a lower-triangular matrix
    % b must be a vector of the same leading dimension as L

    n = size(L, 1);
    y = zeros(n, 1);
    for i = 1:n
        tmp = b(i);
        for j = 1:i-1
            tmp = tmp - L(i,j) * y(j);
        end
        y(i) = tmp / L(i,i);
    end
end
