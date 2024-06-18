C=7e-3;
%Ts=1e-4;
Ts=1e-4;
A=[1 0 0 -1 0 0; 1 0 0 0 -1 0; 1 0 0 0 0 -1; 0 1 0 -1 0 0; 0 1 0 0 -1 0; 0 1 0 0 0 -1; 0 0 1 -1 0 0; 0 0 1 0 -1 0; 0 0 1 0 0 -1;]';
Ecref=(1.0e+07)*3.2781*ones(9,1);
E0=1.0e+06*3.1606*ones(9,1);
lambda=40;
v0 = [-271.0588,-48.6729,294.3140];
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

vchi = reshape(vxyp,[6*Np,1]);
ichi=reshape(ixypred,[9*Np,1]);

Aeq = [];
for i=1:(Np)
    Aeq= cat(2,Aeq,A);
end

beq= zeros(6,1);




[H1,b1]=RealModeliz(E0,Ecref,Ts,lambda,Np,v0,vchi,ichi,A,Q);

N1 = [1 1 1 1; -1 0 -1 0; 0 -1 0 -1; -1 -1 0 0; 1 0 0 0; 0 1 0 0; 0 0 -1 -1; 0 0 1 0; 0 0 0 1];

%Creando la matriz N diagonal por bloques
N = zeros([9, (Np)*4]);

for i=1:Np
    N((i-1)*size(N1,1)+1:i*size(N1,1), (i-1)*size(N1,2)+1:i*size(N1,2)) = N1;
end

%Matrices transformadas
H = (N')*H1*N;
b = (N')*b1;

[~,n2]=size(H1);
[~,n]= size(H);
Ilb=-250*ones(9*Np,1);
Iub=250*ones(9*Np,1);

Bdesg = [-N;N];
bdesg= [-Ilb;Iub];
X0 = InitialCondition(-Bdesg,-bdesg);

%Inf es el mínimo cantidad de iteraciones que queremos que hagan y sup es
%el máximo
inf = 1;
sup = 25;
m= sup - inf + 1;

fvalIPPropio=zeros(m,1);
tiempoIPPropio=zeros(m,1);

tiempoIPLinSolve=zeros(m,1);
fvalIPLinSolve=zeros(m,1);

tiempoIPFminCon=zeros(m,1);
fvalIPNewFminCon=zeros(m,1);



for j= inf:sup
    i=j-(inf-1);
    %Aquí usamos Interior Point con nuestro método
    tic;
    [~, fval] = Interior_Point_Propio(2*H,b,-Bdesg,-bdesg,X0,j);
    tiempoIPPropio(i)=toc;

    fvalIPPropio(i) = fval;

    %Aquí usamos Interior Point con Linsolve
    tic;
    [~, fval2] = IntPointLinSolve(2*H,b,-Bdesg,-bdesg,X0,j);
    tiempoIPLinSolve(i) = toc;
    fvalIPLinSolve(i)=fval2;

    %Aquí usamos Interior Point con Fmincon
    options = optimoptions('fmincon', 'Algorithm', 'interior-point', 'MaxIterations', j, 'Display', 'off');
    tic;
    [~, fval3] = fmincon(@(x)quadratic(x,H,b),X0,Bdesg,bdesg,[],[],[],[],[], options);
    tiempoIPFminCon(i) = toc;
    fvalIPFminCon(i) =fval3;
end



subplot(2, 1, 1); % 2 filas, 1 columna, primer subgráfico

plot(inf:sup, tiempoIPPropio, 'r');

hold on;
plot(inf:sup, tiempoIPLinSolve, 'b');

hold on;
plot(inf:sup,tiempoIPFminCon)

legend('IP','LinSolve', 'Fmincon');

title('Comparación de métodos');
xlabel('Número de iteraciones');
ylabel('Tiempo de ejecución');

figure;

subplot(2, 1, 2); % 2 filas, 1 columna, segundo subgráfico
plot(inf:sup, fvalIPPropio, 'r', 'LineWidth', 2, 'LineStyle', '--');

hold on;
plot(inf:sup, fvalIPLinSolve, 'b');

hold on;
plot(inf:sup, fvalIPFminCon);

legend('IP','LinSolve', 'Fmincon');

title('Comparación de métodos');
xlabel('Número de iteraciones');
ylabel('Valor de las funciones');