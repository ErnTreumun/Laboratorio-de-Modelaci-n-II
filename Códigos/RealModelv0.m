function [H,b,c] =RealModelv0(E0,Eref, Ts, lambda0, Np, iz, vchi,ichi, A,Q,lambdaz)
    [m,n] = size(A); %m=6, n=9;
    %Aquí creamos una función tal que, dado el i_z, las restricciones de caja
    %de v, y también todas las demás matrices de información, nos arroja las
    %matrices A y b para que el problema quede escrito de la forma
    % Minimizar v_0^T H v_0 + b^T v_0 s.t. (box constraints).
    H = zeros(Np, Np);
    b = zeros(Np,1);
    KL = zeros(n,1);
    XLaux = [];
    c = lambdaz*(iz'*iz);
    for L = 1: Np
        xL = ichi(n*(L-1)+1: n*L) + iz(n*(L-1)+1: n*L);
        %PL tiee tamaño 9
        PL = Ts * [diag(A' * vchi(m*(L-1)+1: m*L)) * xL];
        KL = KL + PL;
        ML = E0 + KL - Eref;

        c = c + (ML') * Q * ML;

        XLaux = [XLaux, xL];
        XL = cat(2, XLaux, zeros(n,Np - L));
        
        H = H + XL'*Q*XL;
        b = b + 2 * ( (ML)'* Q * XL)';
    end
    H = H + lambda0 * eye(Np);