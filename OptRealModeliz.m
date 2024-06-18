function [HatQ,b] =OptRealModeliz(E0,Eref, Ts, lambda, Np, v0, vchi,ichi, A,Q)
    [m,n]=size(A); %m=6, n=9;
    %A=[1 0 0 1 0 0 1 0 0; 0 1 0 0 1 0 0 1 0; 0 0 1 0 0 1 0 0 1; -1 -1 -1 0 0 0 0 0 0; 0 0 0 -1 -1 -1 0 0 0; 0 0 0 0 0 0 -1 -1 -1];
    b= zeros(n*(Np),1);
    HatQ = zeros((Np)*n,(Np)*n);
    CL=zeros(n,1);
    HatB = [];
    for L=1: Np
        %Aquí ichi tiene tamaño n*(Np)
        BL = Ts* [diag(A'* vchi(m*(L-1)+1: m*L)) + v0(L)*eye(n)];
        HatB = [HatB,BL];
        CL = CL + BL* ichi(n*(L-1)+1: n*L);
        WL = E0 + CL - Eref;
        b= b + cat(1,2* (HatB' *Q)*WL, zeros(n*(Np)-n*L,1));

		%Cálculo de la matriz Hiz
		for j= L:Np
			%Elemento matriz B(j)
			Bj = Ts* [diag(A'* vchi(m*(j-1)+1: m*j)) + v0(j)*eye(n)];
			%Cálculo de x(i)^\top Q x(j)
			value = Bj'* Q * BL;
			coef = Np - L + 1 - (j-L);
            % Llenar el bloque correspondiente en HatQ
            HatQ((L-1)*n + 1:L*n, (j-1)*n + 1:j*n) = coef * value;
            HatQ((j-1)*n + 1:j*n, (L-1)*n + 1:L*n) = coef * value; % Usando la simetría
		end
    end
    HatQ = HatQ + lambda*eye(n*(Np));