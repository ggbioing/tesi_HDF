function [TT, YY] = emodiafilt(paziente)

global Donnan hind rho Vrc Tp0 Tpis0 Vpl0 Vis0 Pac0 Pis0 alpha Rp gamma ...
       kf sigma k beta Gic type_hdf  eff Td Quf Qf Qb  CL Cd Cs Qs Pvc  ...
       LP La Lv Pn_eq Cc Eis soluto 

eval(['load ',char(paziente)])    

%--------------------------------------------------------------------------
%- INIZILIAZZAZIONE DELLE COSTANTI ----------------------------------------
%--------------------------------------------------------------------------
Tamb    =  25;        % [grad. C]                      Temperatura ambiente
kf      = .004;       % [L^2/mmol/sec]   perm. della membr. cell. all'acqua
gamma   = .93;        % [adim]                     coefficiente di attività
Eis     =  2.45;      % [mmHg/L]       elastanza del compart. interstiziale
Pis0    = -5.97;      % [mmHg]             pressione interstiziale iniziale
Pac0    = 35;         % [mmHg]    pressione iniziale ai capillari arteriosi
Pvc     = 15;         % [mmHg]      pressione ai capillari venosi: costante
rho     = RHO;        % [adim] coef. di permeabilità relativa dei capillari
La      = 3e-4/60;    % [L/mmHg/sec]     coef. di perm. capillari arteriosi
Lv      = 1.86e-3/60; % [L/mmHg/sec]     coef. di perm. capillari    venosi
rHF80s  = 1/(1-sqrt((1-sqrt(.56))))*1470; % [pm]   raggio pori filtro HF80s

%--------------------------------------------------------------------------
%- FORMULE UTILI ----------------------------------------------------------
%--------------------------------------------------------------------------
Rp = @(x) (1-.0107.*x); % [gr/dL]->[adim]      frazione di acqua nel plasma
LP = @(c) 2.1*c + .16*c.^2 + .009*c.^3;    % formula di Landis-Pappenheimer
Staverman = @(r) (1-(1-r/rHF80s).^2).^2;         % calcolo coeff. Staverman

%--------------------------------------------------------------------------
%- DATI E PARAMETRI DEI SOLUTI --------------------------------------------
%--------------------------------------------------------------------------
Tpis0  = .2/1.2*Tp0;        % [gr/L] proteine interstiziali totali iniziali
soluto = {'Na+',   'K+',  'Cl-', 'Ca++', 'P03-', 'Mg++', 'urea',  'creat'};
MW     = [ 23       39      35     40      96     24       60       113  ];
rs     = [ 100     140     180    100     238     72      160       300  ];
alpha  = [139/142 4/4.2 108/108 1.2/1.3  2.3/2  .7/.8     4/4     .2/.2  ]; 
beta   = [ 14/139 140/4   4/108   eps   11/2.3  20/.7     4/4      9/.2  ];
Donnan = @(Tp) (alpha-1)'*Tp/7 +1; %coeff. di Donnan in funz. di alpha e Tp
k      = PAR(:,1) ;       % coeff. trasferimento di massa alla membr. cell.
sigma  = ones(8,1);               % coeff. di riflessione alla membr. cell.
hind   =   1 - Staverman(rs);  % coeff. di hindrance della membr. capillare
eff    = PAR(:,2)  ;         % parametro di ``efficienza'' del dializzatore
Gic    = zeros(8,1);                  % tasso di generazione intracellulare

%--------------------------------------------------------------------------
%- IMPOSTAZIONI INIZIALI MACCHINA DIALIZZATRICE (FRESENIUS 5008S) ---------
%--------------------------------------------------------------------------
Quf = Vuf/Td/60;                      % [L/sec] portata di ultrafiltrazione
Qf  = Qs + Quf ;                      % [L/sec] portata di filtrazione

%- Composizione delle sacche ACIDE contenente Sodio e altri soluti (mmol/L)
%------ 'Na', 'K', 'Cl',  'Ca', 'P', 'Mg', 'urea', 'creat' ----------------
C_a = [ 4725   90  4995   67.5   0   22.5    0        0  ;... % AX03
        4725  135  5040   67.5   0   22.5    0        0 ];    % A161
if     sacca_m == 'AX03', C_a = C_a(1,:);
elseif sacca_m == 'A161', C_a = C_a(2,:); end

%- Composizione della sacca BASICA contenente Bicarbonato di Sodio (mmol/L)
%- FRESENIUS biBag  (soluzione satura di Bicarbonato di sodio) ------------
Csat_bic = interp1([0 20 60],[69 96 165]/84, Tamb ,'linear','extrap')*1000;
%%      'Na',    'K', 'Cl', 'Ca', 'P', 'Mg', 'urea','creat'
C_b = [ Csat_bic  0    0     0     0     0      0     0  ];

%- Composizione del liquido dializzante e/o di sostituzione ---------------
Qtot = Qd + Qs;               % [L/sec]  flusso totale dopo il miscelamento
Q_b   = Qtot * C_Bic/Csat_bic;                     % [L/sec] portata basica
Q_a   = ( Qtot*C_Na - Q_b*Csat_bic )/C_a(1);       % [L/sec] portata  acida

%- Composizione delle sacche ACIDE contenente Sodio e altri soluti (mmol/L)
%------ 'Na', 'K', 'Cl',  'Ca', 'P', 'Mg', 'urea', 'creat'-----------------
C_a = [ 4725   90  4995   67.5   0   22.5    0        0  ;... % AX03
        4725  135  5040   67.5   0   22.5    0        0 ];    % A161
if     sacca == 'AX03', C_a = C_a(1,:);
elseif sacca == 'A161', C_a = C_a(2,:); end

for i = 1:8
    Cd(i) = ( Q_a*C_a(i) + Q_b*C_b(i) )/Qtot; %[mmol/L] liquido dializzante
end
Cs = Cd

%--------------------------------------------------------------------------
%- CALCOLO DEI VALORI DI CLEAREANCE: Filtro HF80s -------------------------
%--------------------------------------------------------------------------

%-      Urea Creat. B12 Inulina 
mw    = [ 60  113  1355  5200]; % peso molecolare
cl200 = [192  180   135   110]; % clearance a 200 mL/min
cl300 = [248  225   155   120]; % clearance a 300 mL/min

CL200 = interp1(log(mw),cl200,log(MW),'linear','extrap'); % interp. log
CL300 = interp1(log(mw),cl300,log(MW),'linear','extrap'); % interp. log

if type_hdf == 'pre '
  CL=interp1([200 300],[CL200' CL300']',(Qb+Qs)*1000*60,'linear','extrap');
else
  CL=interp1([200 300],[CL200' CL300']',Qb*1000*60,'linear','extrap');
end
CL = CL/1000/60; % [mL/min] -> [L/sec]  Clearance dei soluti

%--------------------------------------------------------------------------
%- CALCOLO DEI VALORI INIZIALI --------------------------------------------
%--------------------------------------------------------------------------
Vtot0 = .6*W0;        % Vtot = 60% del peso corporeo [Guyton et al.]
Vpl0  =  3/42*Vtot0 ; % [ibidem]
Vis0  = 11/42*Vtot0 ; % [ibidem]
Vic0  = 28/42*Vtot0 ; % [ibidem]
Vemo0 = Vpl0/(1-Ht0); % Calcolo del volume ematico in funzione de Vpl e Hct
Vrc   = Vemo0 - Vpl0; % Volume eritrocitario: costante per ipotesi
a_d0  = Donnan(Tp0) ; % coeff. di Donnan in funz. della conc. di proteine

Qic_eq = 0;       % flusso di liquido verso il compartimento intracellulare
for s=1:8
    Cis0(s) = a_d0(s) * Cpl0(s);            % Conc. INTERSTIZIALI  iniziali
    Cic0(s) = beta(s) * Cis0(s);            % Conc. INTRACELLULARI iniziali
    Mic0(s) = Cic0(s)*Vic0;                 %     massa intracell. iniziale
    Mex0(s) = Cis0(s)*(Vpl0/a_d0(s) + Vis0);%        massa extrac. iniziale
    Qic_eq  = Qic_eq + gamma * kf * sigma(s) * ( Cic0(s) - Cis0(s) );
end

Cc=inf; % compliance dei capillari

% Pressione netta di filtrazione (idraulica+colloidosmotica) all'equilibrio
Pn_eq =  ((Pac0+Pvc)/2 - Pis0) - ( LP(Tp0) - LP(Tpis0) );

global Y0 %Vettore contenente i dati iniziali in ingresso al sistema di ODE
Y0 = [Mic0(1:8) Mex0(1:8) Vic0 Vis0 Vpl0 Pac0 Pis0];

%--------------------------------------------------------------------------
%- INTEGRAZIONE NUMERICA DEL SISTEMA DI ODE -------------------------------
%--------------------------------------------------------------------------

%options = odeset('OutputFcn',@odeplot);
[TT,YY] = ode15s(@ODE_emodiafilt,[0:60:60*Td],Y0);

%--------------------------------------------------------------------------
%- RISULTATI DELLA SIMULAZIONE --------------------------------------------
%--------------------------------------------------------------------------
plotResults(TT,YY,'b',char(paziente))
