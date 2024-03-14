%% FLUJO DE POTENCIA DE GAUSS-JACOBI y GAUSS-SEIDEL ACELERADO
%% Bases del sistema
Sbase = 1000; %Potencia base (kVA)
Vbase = 23;   %Tension base (kV)
%% Lineas [i j  R[ohm]  X[ohm]] Informacion de lineas del sistema
Lineas  = [1 2 0.5025 0.3025
           2 3 0.4020 0.2510
           3 4 0.3660 0.1864
           2 5 0.3840 0.1965
           5 6 0.8190 0.7050
           2 7 0.2872 0.4088];
%% Nodos [j  V  P[kW] Q[kvar]] Informacion de nodos del sistema
Nodos  = [1  1 0    0    ;
          2  1 1000 600  
          3  1 900  500  
          4  1 2500 1200 
          5  1 1200 950  
          6  1 1050 780  
          7  1 2000 1150];
%% Equivalente en por unidad del sistema
Zbase         = ((1000*Vbase)^2)/(1000*Sbase);
Lineas(:,3:4) = Lineas(:,3:4)/Zbase; 
Nodos(:,3:4)  = Nodos(:,3:4)/Sbase;
%% Calculo de la matriz de admitancia nodal
N = size(Nodos,1); % Numero de nodos
L = size(Lineas,1);% Numero de lineas
Ybus = zeros(N,N); Ybusx = Ybus;
for l = 1:L
    k = Lineas(l,1); m = Lineas(l,2);
    Zkm = Lineas(l,3)+1i*Lineas(l,4);
    Ybus(k,m) = -1/Zkm; Ybus(m,k) = -1/Zkm;
    Ybusx(k,m) = -1/Zkm; Ybusx(m,k) = -1/Zkm;
    Ybus(k,k) = Ybus(k,k) + 1/Zkm;
    Ybus(m,m) = Ybus(m,m) + 1/Zkm;
end
%% Metodo iterativo de Gauss-Jacobi
tmax = 1000; % Numero de iteraciones maximo
epsilon = 1e-10; % Error de convergencia
alpha =  1; % Factor de aceleracion
V0 = Nodos(:,2); % Voltajes iniciales
Vt = V0; % Vector para calculo de voltajes de la siguiente iteracion
Sd = -(Nodos(:,3) + 1i*Nodos(:,4)); % Vector de demandas del sistema
for t = 1:tmax
    for k = 2:N
        Vt(k,1) = (1/Ybus(k,k))*(conj(Sd(k,1)/V0(k,1)) - Ybusx(k,:)*V0(:));
        Vt(k,1) = V0(k,1) + alpha*(Vt(k,1) - V0(k,1));
    end
    if max(abs(abs(Vt) - abs(V0)))<epsilon
        ploss = real(Vt.'*conj(Ybus*Vt));
        break
    else
        V0 = Vt;
    end
end