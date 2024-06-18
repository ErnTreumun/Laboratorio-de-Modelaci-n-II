function [HatQ,b] =RealModeliz(E0,Eref, Ts, lambda, Np, v0, vchi,ichi, A,Q)
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
        HatBaux = cat(2, HatB, zeros(n, n*(Np)- n*L));
        %Este paso se podría mejorar encontrando explícitamente la multiplicación,
		%debido a la cantidad de ceros que existen. Sin
        %embargo ya que el problema no es tan grande, esto no es tan
		%problemático
        HatQ = HatQ+ HatBaux'*Q*HatBaux;
    end
    HatQ = HatQ + lambda*eye(n*(Np));