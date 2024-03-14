%% FLUJO DE POTENCIA DE LA SECANTE EN VARIABLE COMPLEJA
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
Vdo = ones(l,1)*Vg; 
Sd = (Nodos(2:end,3) + 1j*Nodos(2:end,4) );
Ydd = Ybus(2:end,2:end); 
Ydg = Ybus(2:end,1);
%% Metodo de Broyden
e= 1e-10; % Error
tmax = 100; ex = [];
At = Ydd;
F0 = diag(conj(Vdo))*(Ydd*Vdo + Ydg*Vg) + conj(Sd);
for t = 1:tmax
     St1 = -inv(At)*F0;
     Vd = Vdo + St1;
     if max(abs(abs(Vd) - abs(Vdo))) < e
         Vn = [Vg;Vd];
         break
     else
         Vdo = Vd;
         F1 = diag(conj(Vdo))*(Ydd*Vdo + Ydg*Vg) + conj(Sd);
         Yt1 = F1 - F0; F0 = F1;
         At = At + ((Yt1 - At*St1)*transpose(St1))/(transpose(St1)*St1);
     end
end	
%% Perdidas 
ploss  = real(Vn.'*(conj(Ybus)*conj(Vn)))*Sbase; 	