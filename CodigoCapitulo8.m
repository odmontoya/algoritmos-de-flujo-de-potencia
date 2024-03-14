%% FLUJO DE POTENCIA APROXIMACION DEL PRODUCTO
%% Bases del sistema
Sbase = 100;
Vbase = 4.16;
%% Lineas [i  j  R[ohm] X[ohm]]
Lineas =  [1  2  0.0920 0.0477;
           2  3  0.0493 0.0251;
           3  4  0.0366 0.0186;
           4  5  0.0381 0.0194;
           5  6  0.0819 0.0707;
           3  7  0.0187 0.0618;
           7  8  0.0284 0.0144;
           8  9  0.0591 0.0933;
           9  10 0.0804 0.0700;
           10 11 0.0196 0.0650;
           11 12 0.0745 0.0123;
           9  13 0.0203 0.0103;
           13 14 0.0284 0.0144;
           14 15 0.0598 0.0933];
%% Nodos [ i   V  P[kW] Q[kW]] 
Nodos =  [ 1   1  0     0;
           2   1  100   160;
           3   1  190   140;
           4   1  120   80;
           5   1  100   70;
           6   1  90    50;
           7   1  200   100;
           8   1  200   100;
           9   1  210   110;
           10  1  110   60;
           11  1  150   80;
           12  1  170   110;
           13  1  150   90; 
           14  1  120   80;
           15  1  180   60];
%% Cambio a Por Unidad 
Zbase = ((1000*Vbase)^2)/(1000*Sbase);
Lineas(:,3:4) = Lineas(:,3:4)/Zbase; 
Nodos(:,3:4) = Nodos(:,3:4)/Sbase;
%% Formacion de la Ybus
n = size(Nodos,1); %Numero de nodos 
l = size(Lineas,1); %Numero de lineas
Ybus = zeros(n,n); %Matriz de impedancias del sistema
for i = 1:l
    Ni = Lineas(i,1); Nj = Lineas(i,2);
    ZL = Lineas(i,3) + 1j*Lineas(i,4);
    Ybus(Ni,Nj) = -1/ZL; Ybus(Nj,Ni) = -1/ZL;
    Ybus(Ni,Ni) =  Ybus(Ni,Ni) + 1/ZL;
    Ybus(Nj,Nj) =  Ybus(Nj,Nj) + 1/ZL;
end
%% Asignacion de valores iniciales
Vg = Nodos(1,2); 
Vdo = ones(l,1)*Vg;
Sb = (Nodos(2:end,3) + 1j*Nodos(2:end,4) );
Ydd = Ybus(2:end,2:end); 
Ydg = Ybus(2:end,1);
%% Metodo de aproximacion del producto
e = 1e-10; %Error
tmax = 10000;
for t = 1:tmax
    %%Construccion de las matrices A, B y el vector C.
     A=diag((Ydg*Vg)+(Ydd*Vdo));
     B=diag(conj(Vdo))*Ydd;
     C=conj(Sb)-(diag(conj(Vdo))*Ydd*Vdo);
     %%Matrices y vectores reales e imaginarios
     Ar = real(A); Ai = imag(A);
     Br = real(B); Bi = imag(B);
     Cr = real(C); Ci = imag(C);
     %%Matriz compacta y actualizacion de voltajes
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