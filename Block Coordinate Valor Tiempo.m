C=7e-3;
%Ts=1e-4;
Ts=1;
A=[1 0 0 -1 0 0; 1 0 0 0 -1 0; 1 0 0 0 0 -1; 0 1 0 -1 0 0; 0 1 0 0 -1 0; 0 1 0 0 0 -1; 0 0 1 -1 0 0; 0 0 1 0 -1 0; 0 0 1 0 0 -1;]';
Ecref=(1.0e+07)*3.2781*ones(9,1);
Ag=7000;
Am=5000;
E0=1.0e+06*3.1606*ones(9,1);
lambdaz=1;
lambda0 = 40;
Q= eye(9);


Np=3;
I1=10;
omega1= 300;
omega2=110;
V1=300;
V2=110;

vxyp = 1.0e+03 * [
    0.9992,    0.6432,   -0.1322;
   -0.4642,    0.3415,    0.9245;
   -0.5349,   -0.9847,   -0.7923;
    7.0000,    6.9930,    6.9721;
   -3.5000,   -3.7670,   -4.0265;
   -3.5000,   -3.2261,   -2.9457
];

ixypred = [
   -30.6755,  -29.8695,  -18.4745;
     2.4582,  -16.1386,  -33.1005;
    -4.0280,    4.9888,    1.8632;
   -71.5326,  -65.8675,  -49.5417;
   -38.3989,  -52.1366,  -64.1677;
   -44.8851,  -31.0092,  -29.2040;
    42.4269,   47.1479,   62.3045;
    75.5605,   60.8787,   47.6785;
    69.0743,   82.0061,   82.6422
];

Aeq = [];
for i=1:(Np)
    Aeq= cat(2,Aeq,A);
end

beq= zeros(6,1);

N1 = [1 1 1 1; -1 0 -1 0; 0 -1 0 -1; -1 -1 0 0; 1 0 0 0; 0 1 0 0; 0 0 -1 -1; 0 0 1 0; 0 0 0 1];


N=repmat({N1},1,Np);
N=blkdiag(N{:});


%Matrices para el problema v0
Adesgv0 = [-ones(Np); ones(Np)];
v0max = 5* ones(Np,1);
v0min = -5 * ones(Np,1);

bdesgv0 = [- v0max; v0min];

v0 = zeros(Np,1);

[~,n] = size(A);
n2 = n * Np;
imin=-250*ones(9*Np,1);
imax=250*ones(9*Np,1);

Bdesg = [-N1;N1];
bdesg= [-imin;imax];

%Maximo de iteraciones
maxiterv0 = 12;
maxiteriz = 12;
maxitertotal = 100;

%Condiciones iniciales

%Matrices de desigualdad asociados al problema transformado
Bdesg = [-N;N];
bdesg= [-imin;imax];

v0 = rand(Np,1);
iz0 = InitialCondition(-Bdesg,-bdesg);

for j= 1: maxitertotal
    v0aux = v0;
	iz = iz0;

    tic;
    for i=1:j
	    %Valores para ir actualizando en los métodos

		[Higiz, biz] = RealModeliz(E0, Ecref, Ts, lambdaz, Np, v0aux, reshape(vxyp,[6*Np,1]),reshape(ixypred,[9*Np,1]), A, Q);
		    
	    %Matrices del problema transformado
	    HtransIz = (N')*Higiz*N;
		btransIz = (N')*biz;

	    [izaux, fvaliz] = Interior_Point_Propio(2*HtransIz,btransIz,-Bdesg,-bdesg,iz,maxiteriz);
	    iz = izaux;
		
		[Higv0, bv0, c] = RealModelv0(E0, Ecref, Ts, lambda0, Np, N*iz, reshape(vxyp,[6*Np,1]), reshape(ixypred,[9*Np,1]), A,Q,lambdaz);

		[vnew, fvalv0] = Interior_Point_Propio(2*Higv0, bv0, -Adesgv0, -bdesgv0, v0aux, maxiterv0);
		v0aux = vnew;
    end
    timetotalsinOpt(j) = toc;
	valv0sinOpt(j) = fvalv0 + c;
end

figure;
plot(1:maxitertotal, timetotalsinOpt, '-b');

%title('Comparación de métodos');
xlabel('Número de iteraciones');
ylabel('Tiempo de ejecución');

figure;

plot(1:maxitertotal,valv0sinOpt,  '-b');

xlabel('Número de iteraciones');
ylabel('Valor de la función');