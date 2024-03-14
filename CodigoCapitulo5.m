%% FLUJO DE POTENCIA TRIANGULAR 
clc; clear
%% Bases del sistema
Sbase = 1000;%Potencia base (kVA)
Vbase = 23;%Tension base (kV)
%% Lineas [i j  R[ohm]  X[ohm]] Informacion de lineas del sistema
Lineas =  [1 2  0.5025  0.3025;
           2 3  0.4020  0.2510; 
           3 4  0.3660  0.1864;
           2 5  0.3840  0.1965;
           5 6  0.8190  0.7050;
           2 7  0.2872  0.4088];
%% Nodos [j  V  P[kW] Q[kvar]] Informacion de nodos del sistema
Nodos  = [1   1   0   0;
          2   1  1000 600;
          3   1  900  500;
          4   1  2500 1200;
          5   1  1200 950;
          6   1  1050 780;
          7   1  2000 1150];
%% Equivalente en por unidad del sistema
Zbase         = ((1000*Vbase)^2)/(1000*Sbase);
Lineas(:,3:4) = Lineas(:,3:4)/Zbase; 
Nodos(:,3:4)  = Nodos(:,3:4)/Sbase;
%% Formacion de las matrices  T y Z 
l = size(Lineas,1); % Numero de lineas 
T = zeros(l,l); % Matriz triangular 
Z = zeros(l,l); % Matriz primitiva de impedancias 
Lineas = sortrows(Lineas,2);
for i = l:-1:1
    Z(i,i) = Lineas(i,3) + 1j*Lineas(i,4);
    T(i,i) = 1;
    Ni = Lineas(i,1);
    while Ni ~= 1
        Ni = Ni-1;
        T(Ni,i) = 1;
        Ni = Lineas(Ni,1);   
    end
end
%% Asignacion de valores iniciales 
Vg = Nodos(1,2);
Vd0 = ones(l,1).* Vg;
Zbus = (T.')*Z*T; % Matriz Zbus se crea por facilidad de implementacion
Sd = Nodos(2:end,3) + 1j*Nodos(2:end,4); 
%% Formulacion del metodo triangular
e = 1e-10; 
tmax = 100;
for t=1:tmax
     Vd = ones(l,1)*Vg - Zbus*((inv(diag(conj(Vd0))))*conj(Sd));    
     if max(abs((abs(Vd)-abs(Vd0))))<=e
         Vn = [Vg;Vd];
         break
     else 
         Vd0 = Vd;        
     end 
end
%% Calculo de perdidas 
Id = diag(conj(Vd))\conj(Sd); 
Il = T*Id;
El = Z*Il;
ploss = real(El.'*conj(Il))*Sbase;