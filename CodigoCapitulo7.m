%% FLUJO DE POTENCIA DE APROXIMACION HIPERBOLICA
%% Bases del sistema
Sbase = 1000; % Potencia base (kVA)
Vbase = 13.8; % Tension base (kV)
%% Lineas [i  j    R[ohm] X[ohm]]
Lineas =  [1   2   0.15208  0.19855
           2   3   0.65805  0.59745
           3   4   0.19742  0.17924
           4   5   0.43848  0.26038
           5   6   0.48720  0.28931
           6   7   0.48197  0.22732
           7   8   0.87630  0.41330
           8   9   1.09540  0.51663
           9  10   0.87630  0.41330
           2  11   0.87630  0.41330
           11 12   1.07780  0.50836
           12 13   0.65722  0.30998
           13 14   0.49073  0.23145
           14 15   0.87630  0.41330
           15 16   0.87630  0.41330
           3  17   0.87630  0.41330
           17 18   0.52578  0.24798
           18 19   0.78867  0.37197
           19 20   0.83248  0.39263
           20 21   0.87630  0.41330
           4  22   0.87630  0.41330
           5  23   0.87630  0.41330
           6  24   0.35052  0.16532
           8  25   0.52578  0.24798
           8  26   0.52578  0.24798
           26 27   0.70104  0.33064];
%% Nodos  [j  V  P[kW] Q[kvar]]
Nodos  =  [ 1 1       0        0
            2 1       0        0
            3 1       0        0
            4 1   297.5    184.4
            5 1       0        0
            6 1     255      158
            7 1       0        0
            8 1   212.5    131.7
            9 1       0        0
           10 1   266.1    164.9
           11 1      85     52.7
           12 1     340    210.7
           13 1   297.5    184.4
           14 1   191.3    118.5
           15 1   106.3     65.8
           16 1     255      158
           17 1     255      158
           18 1   127.5       79
           19 1   297.5    184.4
           20 1     340    210.7
           21 1      85     52.7
           22 1   106.3     65.8
           23 1    55.3     34.2
           24 1    69.7     43.2
           25 1     255      158
           26 1    63.8     39.5
           27 1     170    105.4];
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
%% Metodo de aproximacion hiperbolica
e= 1e-10; % Error
tmax = 100;
for t = 1: tmax
    % Construccion de las matrices A, B y el vector C.
    A = diag((1./conj(Vdo)).^2)*diag(conj(Sd));
    B = - Ydd;
    C = - Ydg*Vg - 2*diag(1./conj(Vdo))*conj(Sd);
    % Matrices y vectores reales
    Ar = real(A); Ai = imag(A);
    Br = real(B); Bi = imag(B);
    Cr = real(C); Ci = imag(C);
    % Matriz compacta y actualicion de voltajes
    Jaco = [Ar+Br Ai-Bi;Ai+Bi Br-Ar];
    DeltaC = [Cr;Ci];
    Vri = -Jaco\DeltaC;
    Vd = Vri(1:n-1,1) + 1i*Vri(n:end,1);
     if max(abs(abs(Vd) - abs(Vdo))) < e
         Vn = [Vg;Vd];
         break
     else
         Vdo = Vd;
     end
end	
%% Perdidas 
ploss  = real(Vn.'*(conj(Ybus)*conj(Vn)))*Sbase; 	