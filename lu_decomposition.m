function [P, L, U] = lu_decomposition(A)            %% 8n^3/3-n^2+5n/3 flops
    % Descomposicion PA=LU

    n = length(A);

    % Crear L,U
    L = zeros(n, n);
    U = zeros(n, n);

    P = pivot_matrix(A);                  %% n^2+3n flops
    PA = P * A;                           %% 2n^3-n^2 flops

    % Perform the LU Decomposition
    for j = 1:n                                         %% n^3/3-n^2-4n/3
        % All diagonal entries of L are set to unity
        L(j, j) = 1.0;                       
  
        for i = 1:j                                     %% j^2-j flops
            s1 = sum(U(1:i-1, j) .* L(i, 1:i-1)');      %% 2i-3 flops
            U(i, j) = PA(i, j) - s1;                    %% 1 flop
        end

        for i = j:n                                     %% (2j-1)(n-j+1) flops
            s2 = sum(U(1:j-1, j) .* L(i, 1:j-1)');      %% 2j-3 flops
            L(i, j) = (PA(i, j) - s2) / U(j, j);        %% 2 flops
        end
    end
end

function P = pivot_matrix(M)                                             %% m^2+3m flops
    % Matriz pivote del metodo doolittle
    m = length(M);
    id_mat = eye(m);


    for j = 1:m                             %% m^2+3m flops
        [~, row] = max(abs(M(j:end, j)));   %% 2(m-j+1) flops
        row = row + j - 1;                  %% 2 flops
        if j ~= row
            id_mat([j, row], :) = id_mat([row, j], :);
        end
    end
    P = id_mat;
end
