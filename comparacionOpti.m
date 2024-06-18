C=7e-3;
%Ts=1e-4;
Ts=1;
A=[1 0 0 -1 0 0; 1 0 0 0 -1 0; 1 0 0 0 0 -1; 0 1 0 -1 0 0; 0 1 0 0 -1 0; 0 1 0 0 0 -1; 0 0 1 -1 0 0; 0 0 1 0 -1 0; 0 0 1 0 0 -1;]';
Ecref=C*500^2*ones(9,1)/2;
E0=zeros(9,1);
lambdaz=0.3;
lambda0 = 0.1;
Q= eye(9);


Np=10;
I1=10;
omega1= 300;
omega2=110;
V1=300;
V2=110;

vchi = zeros(6*(Np),1);
ichi=zeros(6*(Np),1);

t=0;
for i = 1:Np
    vchi((i-1)*6+1:i*6,1) = [V1*cos(omega1*t); V1*cos(omega1*t - 2*pi/3); V1*cos(omega1*t + 2*pi/3); V2*cos(omega2*t); V2*cos(omega2*t - 2*pi/3); V2*cos(omega2*t + 2*pi/3)];
    v0((i-1)*9+1:i*9,1) = [0;0;-1;1;2;3;0;2;3];
    ichi=[I1*cos(omega1*t); I1*cos(omega1*t - 2*pi/3); I1*cos(omega1*t + 2*pi/3); (I1/V1)*V2*cos(omega2*t); (I1/V1)*V2*cos(omega2*t - 2*pi/3); (I1/V1)*V2*cos(omega2*t + 2*pi/3)];
    t= t + i*Ts;
end

ichi=ones(9*(Np),1);

Aeq = [];
for i=1:(Np)
    Aeq= cat(2,Aeq,A);
end

beq= zeros(6,1);

N1 = [1 1 1 1; -1 0 -1 0; 0 -1 0 -1; -1 -1 0 0; 1 0 0 0; 0 1 0 0; 0 0 -1 -1; 0 0 1 0; 0 0 0 1];

%Creando la matriz N diagonal por bloques
N = zeros([9, (Np)*4]);

for i=1:Np
    N((i-1)*size(N1,1)+1:i*size(N1,1), (i-1)*size(N1,2)+1:i*size(N1,2)) = N1;
end



%Matrices para el problema v0
Adesgv0 = [-zeros(Np); zeros(Np)];
v0max = 5* ones(Np,1);
v0min = -5 * ones(Np,1);

bdesgv0 = [- v0max; v0min];

v0 = zeros(Np,1);


[~,n] = size(A);
n2 = n * Np;
imin=-5*ones(n2,1);
imax=5*ones(n2,1);


%Maximo de iteraciones
maxiterv0 = 10;
maxiteriz = 10;
maxitertotal = 30;



%Vamos a usar esto para analizar el código
valiz = zeros(maxitertotal,1);
valv0 = zeros(maxitertotal,1);
timetotal = zeros(maxitertotal, 1);

%Condiciones iniciales
v0 = zeros(Np,1);
iz0 = zeros(4*Np,1);

%Ahora sin haber optimizado
valv0sinOpt = zeros(maxitertotal,1);
timetotalsinOpt= zeros(maxitertotal, 1);

for j= 1: maxitertotal
	%Valores para ir actualizando en los métodos
	v0aux = v0;
	iz = iz0;

	tic;
	for i=1:j
		%Calcula la matriz para el problema no transformado
		[Higiz, biz] = RealModeliz(E0, Ecref, Ts, lambdaz, Np, v0aux, vchi, ichi, A, Q);
		
		%Matrices del problema transformado
		HtransIz = (N')*Higiz*N;
		btransIz = (N')*biz;
		
		%Matrices de desigualdad asociados al problema transformado
		Bdesg = [-N;N];
		bdesg= [-imin;imax];

		[izaux, fvaliz] = IntPointf(2*HtransIz,btransIz,-Bdesg,-bdesg,iz,maxiteriz);
		iz = izaux;
		
		[Higv0, bv0, c] = RealModelv0(E0, Ecref, Ts, lambda0, Np, N*iz, vchi, ichi, A,Q,lambdaz);

		[vnew, fvalv0] = IntPointf(2*Higv0, bv0, Adesgv0, bdesgv0, v0aux, maxiterv0);
		v0aux = vnew;
	end
	timetotalsinOpt(j) = toc;
	valv0sinOpt(j) = fvalv0 + c;
	
end

for j= 1: maxitertotal
	%Valores para ir actualizando en los métodos
	v0aux = v0;
	iz = iz0;

	tic;
	for i=1:j
		%Calcula la matriz para el problema no transformado
		[Higiz, biz] = OptRealModeliz(E0, Ecref, Ts, lambdaz, Np, v0aux, vchi, ichi, A, Q);
		
		%Matrices del problema transformado
		HtransIz = (N')*Higiz*N;
		btransIz = (N')*biz;
		
		%Matrices de desigualdad asociados al problema transformado
		Bdesg = [-N;N];
		bdesg= [-imin;imax];

		[izaux, fvaliz] = IntPointf(2*HtransIz,btransIz,-Bdesg,-bdesg,iz,maxiteriz);
		iz = izaux;
		
		[Higv0, bv0, c] = OptRealModelv0(E0, Ecref, Ts, lambda0, Np, N*iz, vchi, ichi, A,Q,lambdaz);

		[vnew, fvalv0] = IntPointf(2*Higv0, bv0, Adesgv0, bdesgv0, v0aux, maxiterv0);
		v0aux = vnew;
	end
	timetotal(j) = toc;
	valv0(j) = fvalv0 + c;
	
end

figure;
plot(1:maxitertotal, timetotal, '-r', 'DisplayName', 'Optimizado');
hold on;
plot(1:maxitertotal, timetotalsinOpt, '-b', 'DisplayName', 'No Optimizado');
hold off;

%title('Comparación de métodos');
xlabel('Número de iteraciones');
ylabel('Tiempo de ejecución');

legend;

figure;

plot(1:maxitertotal,valv0,  '--r','LineWidth', 2, 'DisplayName', 'Optimizado');
hold on;
plot(1:maxitertotal,valv0sinOpt,  '-b', 'DisplayName', 'No Optimizado');
hold off;

xlabel('Número de iteraciones');
ylabel('Valor de la función');

legend;

