clc
clear

%% Datos de la matriz y las fibras

%Fibra de vidrio (Isótropo)
D_f1=2.6;               %g/cm3
E_f1=72;                %GPa
v_f1=0.2;
G_f1=E_f1/(2*(1+v_f1)); %GPa

%Fibra de carbono
D_f2=1.58;              %g/cm3
E_f21=220;              %GPa
E_f22=22;               %GPa
G_f2=25;                %GPa
v_f2=0.15;

%Resina epóxica (Isótropo)
D_m=1.1;                %g/cm3
E_m=3.3;                %GPa
v_m=0.35;
G_m=E_m/(2*(1+v_m));    %GPa

V_f=0.5;                %Fracción de fibra
D_p1=220;               %Densidad planar [g/m^3]
D_p2=440;
N=8;                    %Número de capas

%% Propiedades material 1 (Fibra de vidrio)

h_11=((D_p1/1e+4)*6)/(D_f1*V_f)*10        %Espesor [mm]
h_12=((D_p2/1e+4)*2)/(D_f1*V_f)*10        %Espesor [mm]

E_11=E_f1*V_f+E_m*(1-V_f)               %Módulo elástico longitudinal [GPa]
E_12=E_m*(1/((1-V_f)+E_m/E_f1*V_f))     %Módulo elástico transversal [GPa]

v_112=v_f1*V_f+v_m*(1-V_f)              %Coeficiente de Poisson 12
v_121=v_112*E_12/E_11                   %Coeficiente de Poisson 21
G_112=G_m*(1/((1-V_f)+G_m/G_f1*V_f))    %Módulo de corte [GPa]

%% Propiedades material 2 (Fibra de carbono)

h_21=((D_p1/1e+4)*6)/(D_f2*V_f)*10        %Espesor [mm]
h_22=((D_p2/1e+4)*2)/(D_f2*V_f)*10        %Espesor [mm]

E_21=E_f21*V_f+E_m*(1-V_f)              %Módulo elástico longitudinal [GPa]
E_22=E_m*(1/((1-V_f)+E_m/E_f22*V_f))    %Módulo elástico transversal [GPa]

v_212=v_f2*V_f+v_m*(1-V_f)              %Coeficiente de Poisson 12
v_221=v_212*E_22/E_21                   %Coeficiente de Poisson 21
G_212=G_m*(1/((1-V_f)+G_m/G_f2*V_f))    %Módulo de corte [GPa]

%% Cálculo de la matriz Q (material 1)

Q_1=[E_11/(1-v_112*v_121)       v_112*E_12/(1-v_112*v_121)  0;
     v_112*E_12/(1-v_112*v_121) E_12/(1-v_112*v_121)        0;
     0                          0                           G_112]

%% Cálculo de la matriz Q (material 2)

Q_2=[E_21/(1-v_212*v_221)       v_212*E_22/(1-v_212*v_221)  0;
     v_212*E_22/(1-v_212*v_221) E_22/(1-v_212*v_221)        0;
     0                          0                           G_212]

%% Matrices T Y R

theta=[0 -pi/4 +pi/4 pi/2];

for i=1:4
m(i)=cos(theta(i));
n(i)=sin(theta(i));

T{i}=[m(i)^2        n(i)^2      2*m(i)*n(i);
      n(i)^2        m(i)^2      -2*m(i)*n(i);
      -m(i)*n(i)    m(i)*n(i)   m(i)^2-n(i)^2];
end

R=[1 0 0; 0 1 0; 0 0 1/2];

%% Matriz de rigidez (material 1)

for i=1:4
Q{1,i}=T{1,i}^-1*Q_1*R*T{1,i}*R^-1;
end

A_1=2*h_12*Q{1,1}+2*h_11*(Q{1,2}+Q{1,3}+Q{1,4})

%% Matriz de rigidez (material 2)

for i=1:4
Q{2,i}=T{1,i}^-1*Q_2*R*T{1,i}*R^-1;
end

A_2=2*h_22*Q{2,1}+2*h_21*(Q{2,2}+Q{2,3}+Q{2,4})

%% Propiedades efectivas

E_a1=1/(6*h_11+2*h_12)*(A_1(1,1)-A_1(1,2)^2/A_1(2,2))
v_a1=A_1(1,2)/A_1(1,1)
G_a1=1/(6*h_11+2*h_12)*A_1(3,3)

E_a2=1/(6*h_21+2*h_22)*(A_2(1,1)-A_2(1,2)^2/A_2(2,2))
v_a2=A_2(1,2)/A_2(1,1)
G_a2=1/(6*h_21+2*h_22)*A_2(3,3)
