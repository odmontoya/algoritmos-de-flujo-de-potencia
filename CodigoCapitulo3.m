%% FLUJO DE POTENCIA BARRIDO ITERATIVO 
%% Bases del sistema
Sbase = 1000;%Potencia base (kVA)
Vbase = 23;%Tension base (kV)
%% Lineas [i j  R[ohm]  X[ohm]] Informacion de lineas del sistema
Lineas =  [1  2 0.3967 0.5290;  
           2  3 0.4232 0.5819;  
           2  4 0.4761 0.9522;  
           4  5 0.2116 0.2116;  
           1  6 0.5819 0.5819;  
           6  7 0.4232 0.5819;  
           6  8 0.5819 0.5819;  
           7  9 0.5819 0.5819;  
           7  10 0.4232 0.5819; 
           1  11 0.5819 0.5819; 
           11 12 0.4761 0.6348; 
           11 13 0.4232 0.5819; 
           13 14 0.2116 0.2116];
%% Nodos [j  V  P[kW] Q[kvar]] Informacion de nodos del sistema
Nodos =  [1  1 0    0    ;
          2  1 2000 1600 ; 
          3  1 3000 400  ; 
          4  1 2000 -400 ; 
          5  1 1500 1200 ; 
          6  1 4000 2700 ; 
          7  1 5000 1800 ; 
          8  1 1000 900  ; 
          9  1 600  -400  ; 
          10 1 4500 -1700; 
          11 1 1000 900  ; 
          12 1 1000 -1100; 
          13 1 1000 900  ; 
          14 1 2100 -800];
%% Equivalente en por unidad del sistema
Zbase         = ((1000*Vbase)^2)/(1000*Sbase);
Lineas(:,3:4) = Lineas(:,3:4)/Zbase; 
Nodos(:,3:4)  = Nodos(:,3:4)/Sbase;
%% Formacion de las matrices A y Z
n = size(Nodos,1);%Numero de nodos
l = size(Lineas,1);%Numero de lineas
A = zeros(n,l);
Z = zeros(l,l);%Matriz primitiva de impedancias
for i=1:l
    Ni = Lineas(i,1); Nj=Lineas(i,2);
    A(Ni,i) =  1;
    A(Nj,i) = -1;
    Z(i,i) = Lineas(i,3)+1j*Lineas(i,4);%Formacion de la matriz primitiva de impedancias
end
%% Asignacion de valores iniciales
Ag = A(1,:); Ad = A(2:end,:); %Separacion de la matriz A
Ydd = Ad*inv(Z)*Ad.'; Ydg = Ad*inv(Z)*Ag.';
Vdo = Nodos(2:end,2);%Tensiones iniciales de los nodos
Vg = Nodos(1,2);%Tension del nodo slack
Sd = Nodos(2:end,3) + 1j*Nodos(2:end,4);
%% Formulacion del metodo de barrido iterativo
e=1e-10;
tmax=100;
for t = 1:tmax
    Vd = -inv(Ydd)*((inv(diag(conj(Vdo))))*conj(Sd)+Ydg*Vg);
    if max(abs((abs(Vd)-abs(Vdo))))<=e
        V=[Vg;Vd];
        break
    else
        Vdo=Vd;
    end
end
%% Calculo de perdidas
Er=Ag.'*Vg+Ad.'*Vd;
Jr= Z\Er;
ploss=real(sum((Z*(abs(Jr)).^2)));