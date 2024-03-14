%% FLUJO DE POTENCIA APROXIMACIONES SUCESIVAS
%% Bases del sistema
Sbase = 10; % Potencia base (kVA)m
Vbase = 23; % Tension base (kV)
%% Lineas [i  j    R[ohm] X[ohm]]
Lineas =  [1  2   0.1233  0.4127; 2  3   0.0140  0.6051; 3  4   0.7463  1.2050; 
           4  5   0.6984  0.6084; 5  6   1.9831  1.7276; 6  7   0.9053  0.7886;
           7  8   2.0552  1.1640; 8  9   4.7953  2.7160; 9  10  5.3434  3.0260];
%% Nodos  [j  V  P[kW] Q[kvar]]
Nodos  =  [1  1  0     0;    2  1  1840  460;
           3  1  980   340;  4  1  1790  446;
           5  1  1598  1840; 6  1  1610  600;
           7  1  780   110;  8  1  1150  60;
           9  1  980   130;  10 1  1640  200];
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
Vdo = ones(l,1)*Vg; 
Sd = (Nodos(2:end,3) + 1j*Nodos(2:end,4) );
Ydd = Ybus(2:end,2:end); 
Ydg = Ybus(2:end,1);
%% Metodo de aproximaciones sucesivas
e= 1e-10; % Error
tmax = 100;
for t = 1: tmax
     Vd = -inv(Ydd)*(inv(diag(conj(Vdo)))*conj(Sd)+ Ydg*Vg);
     if max(abs(abs(Vd) - abs(Vdo))) < e
         Vn = [Vg;Vd];
         break
     else
         Vdo = Vd;
     end
end	
%% Perdidas 
ploss  = real(Vn.'*(conj(Ybus)*conj(Vn)))*Sbase; 	