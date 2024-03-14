%% FLUJO DE POTENCIA DE LA SECANTE EN VARIABLE REAL
%% Bases del sistema
Sbase = 100; % Potencia base (kVA)
Vbase = 4.16; % Tension base (kV)
%% Lineas [i  j    R[ohm] X[ohm]]
Lineas =  [1  2   0.0123  0.0412; 
           2  3   0.0140  0.0605;
           2  4   0.0746  0.2050;
           4  5   0.0698  0.0608];
%% Nodos  [j  V  P[kW] Q[kvar]]
Nodos  =  [1  1  0     0;
           2  1  1750 1300;
           3  1  1525 1400;
           4  1  1375 1250;
           5  1  1400 1200;
           ];
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
%% Metodo de Broyden
e= 1e-10; % Error
tmax = 100; ex = [];
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
Jaco = [H M;N O]; At = Jaco;
Pcal = zeros(n-1,1); Qcal = zeros(n-1,1);
for i = 2:n
    for j = 1:n
        Pcal(i-1,1) = Pcal(i-1,1) + abs(Ybus(i,j))*V0(i,1)*V0(j,1)*cos(a0(i,1) - a0(j,1) - angle(Ybus(i,j)));
        Qcal(i-1,1) = Qcal(i-1,1) + abs(Ybus(i,j))*V0(i,1)*V0(j,1)*sin(a0(i,1) - a0(j,1) - angle(Ybus(i,j)));
    end
end
Dp = Pcal + Pd;  Dq = Qcal + Qd;
Dpq = [Dp;Dq];
F0 = Dpq;
for t = 1: tmax
    St1 = -inv(At)*F0;
    x = [a0(2:n,1);V0(2:n,1)] + St1;
    Vd = V0; a = a0;
    a(2:n,1) = x(1:n-1,1); Vd(2:n,1) = x(n:end,1);
    if max(abs(abs(Vd) - abs(V0))) < e
        Vn  = Vd.*cos(a) + 1i*Vd.*sin(a);
        break
    else
        V0 = Vd; a0 = a;
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
        F1 = Dpq;
        Yt1 = F1 - F0; F0 = F1;
        At = At + ((Yt1 - At*St1)*transpose(St1))/(transpose(St1)*St1);
    end
end
%% Perdidas
ploss  = real(Vn.'*(conj(Ybus)*conj(Vn)))*Sbase;