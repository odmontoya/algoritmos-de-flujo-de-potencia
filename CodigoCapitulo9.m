%% FLUJO DE POTENCIA DE CUASI-NEWTON
%% Bases del sistema
Sbase = 100; % Potencia base (kVA)
Vbase = 10; % Tension base (kV)
%% Lineas [i  j    R[ohm] X[ohm]]
Lineas =  [1  2   0.0122  0.0286; 
           2  3   0.0140  0.0605;
           3  4   0.0355  0.0482;
           4  5   0.0698  0.0608;
           5  6   0.0656  0.0804;
           6  7   0.0901  0.1500;   
           7  8   0.0756  0.0950;
           8  9   0.0127  0.0388;
           9  10  0.0545  0.0836;
           3  11  0.0165  0.0308;
           11 12  0.0185  0.0544;
           12 13  0.0854  0.0102;
           7  14  0.0662  0.0356;
           14 15  0.0689  0.0756;
           9  16  0.0788  0.0968;
           16 17  0.0522  0.0756;
           5  18  0.0165  0.0408;
           18 19  0.0185  0.0644;
           19 20  0.0854  0.0202;
           10 21  0.0662  0.0256; 
           21 22  0.0689  0.0756;
           22 23  0.0788  0.0868];
%% Nodos  [j  V  P[kW] Q[kvar]]
Nodos  =  [1  1  0     0;
           2  1  500   300;
           3  1  700   450;
           4  1  850   600;
           5  1  450   200;
           6  1  800   600;
           7  1  750   300;
           8  1  1000  750;
           9  1  850   500;
           10 1  960  -500;
           11 1  600   400;
           12 1  450   250;
           13 1  400  -800;
           14 1  150   300;
           15 1  725  -1000;
           16 1  680   420;
           17 1  720   560;
           18 1  600   400;
           19 1  450   250;
           20 1  800   400;
           21 1  150   300;
           22 1 1000   725;
           23 1  350   150];
%% Equivalente en por unidad del sistema
Zbase         = ((1000*Vbase)^2)/(1000*Sbase);
Lineas(:,3:4) = Lineas(:,3:4)/Zbase;
Nodos(:,3:4)  = Nodos(:,3:4)/Sbase;
%% Formacion de la Ybus
n = size(Nodos,1); % Numero de nodos
l = size(Lineas,1); % Numero de lineas
Ybus= zeros(n,n); % Matriz de admitancias del sistema
for i = 1:l
    Ni = Lineas(i,1); Nj = Lineas(i,2);
    ZL = Lineas(i,3) + 1j*Lineas(i,4);
    Ybus(Ni,Nj) = -1/ZL;
    Ybus(Nj,Ni) = -1/ZL;
    Ybus(Ni,Ni) =  Ybus(Ni,Ni) + 1/ZL;
    Ybus(Nj,Nj) =  Ybus(Nj,Nj) + 1/ZL;
end
%% Asignacion de valores iniciales
Vg = Nodos(1,2);
V0 = ones(n,1)*Vg; a0 = zeros(n,1);
Pd = Nodos(2:end,3); Qd = Nodos(2:end,4);
%% Metodo de cuasi-Raphson
e= 1e-10; % Error
tmax = 100;
% Construccion del Jacobiano (matrices H, M, N y O)
H = zeros(n,n); M = zeros(n,n); N = zeros(n,n); O = zeros(n,n);
for i = 1:n
    for j = 1:n
        if j ~= i
            % Matriz H
            H(i,j) = abs(Ybus(i,j))*V0(i,1)*V0(j,1)*sin(a0(i,1) - a0(j,1) - angle(Ybus(i,j)));
            H(i,i) = H(i,i) - abs(Ybus(i,j))*V0(i,1)*V0(j,1)*sin(a0(i,1) - a0(j,1) - angle(Ybus(i,j)));
            % Matriz M
            M(i,j) = abs(Ybus(i,j))*V0(i,1)*cos(a0(i,1) - a0(j,1) - angle(Ybus(i,j)));
            M(i,i) = M(i,i) + abs(Ybus(i,j))*V0(j,1)*cos(a0(i,1) - a0(j,1) - angle(Ybus(i,j)));
            % Matriz N
            N(i,j) = -abs(Ybus(i,j))*V0(i,1)*V0(j,1)*cos(a0(i,1) - a0(j,1) - angle(Ybus(i,j)));
            N(i,i) = N(i,i) + abs(Ybus(i,j))*V0(i,1)*V0(j,1)*cos(a0(i,1) - a0(j,1) - angle(Ybus(i,j)));
            % Matriz O
            O(i,j) = abs(Ybus(i,j))*V0(i,1)*sin(a0(i,1) - a0(j,1) - angle(Ybus(i,j)));
            O(i,i) = O(i,i) + abs(Ybus(i,j))*V0(j,1)*sin(a0(i,1) - a0(j,1) - angle(Ybus(i,j)));
        end
    end
    M(i,i) = M(i,i) + 2*real(Ybus(i,i))*V0(i,1);
    O(i,i) = O(i,i) - 2*imag(Ybus(i,i))*V0(i,1);
end
Borrar = 1; % Se elimina la informacion del nodo fuente
H(:,Borrar) = []; H(Borrar,:) = [];
M(:,Borrar) = []; M(Borrar,:) = [];
N(:,Borrar) = []; N(Borrar,:) = [];
O(:,Borrar) = []; O(Borrar,:) = [];
Jaco = [H M;N O]; JacoInv = inv(Jaco);
for t = 1: tmax
    % Potencias calculadas
    Pcal = zeros(n-1,1); Qcal = zeros(n-1,1);
    for i = 2:n
        for j = 1:n
            Pcal(i-1,1) = Pcal(i-1,1) + abs(Ybus(i,j))*V0(i,1)*V0(j,1)*cos(a0(i,1) - a0(j,1) - angle(Ybus(i,j)));
            Qcal(i-1,1) = Qcal(i-1,1) + abs(Ybus(i,j))*V0(i,1)*V0(j,1)*sin(a0(i,1) - a0(j,1) - angle(Ybus(i,j)));
        end
    end
    Dp = Pcal + Pd;  Dq = Qcal + Qd;
    % Se calcula el jacobiano reducido
    Dpq = [Dp;Dq];
    Dx = JacoInv*Dpq;
    V = V0;
    V(2:n,1) = V0(2:n,1) - Dx(n:end,1);
    a = a0;
    a(2:n,1) = a0(2:n,1) - Dx(1:n-1,1);
    if max(abs(V - V0)) < e
        Vn = V.*cos(a) + 1i*V.*sin(a);
        break
    else
        V0 = V;    a0 = a;
    end
end
%% Perdidas
ploss  = real(Vn.'*(conj(Ybus)*conj(Vn)))*Sbase;