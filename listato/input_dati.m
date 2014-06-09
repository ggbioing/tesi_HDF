nome    = '**********';
cogmone = '*****'     ;
ID      = 'AMB1'      ; % identificativo paziente\seduta

%--------------------------------------------------------------------------
%-------------------- PRELIEVI EMATICI ------------------------------------
%--------------------------------------------------------------------------
t_sample = [      0       60     120     180     240    ]'; % min
 
urea  =    [     84       50     35      26      20     ];  % mg/dL
gluc  =    [     141       91    74      72      74     ];  % mg/dL
creat =    [    7.46     4.26   3.29    2.64    2.31    ];  % mg/dL
Na    =    [     142      140    140     140     140    ];      % mEq/L
K     =    [     4.6       4      4      3.8     3.8    ];      % mEq/L
Cl    =    [     101      104    103     103     102    ];      % mEq/L
Ca    =    [     8.6      9      9.2     9.6     9.6    ];  % mg/dL
P     =    [     4.7      2.1    1.9     1.9     1.9    ];  % mg/dL
Mg    =    [     1.8      1.8    1.8     1.7     1.7    ];  % mg/dL
Tp    =    [     6.6      6.6    6.8     7.1     7.3    ];      % gr/dL
alb   =    [    3432     3419   3592    3666    3903    ];  % mg/dL
%--------------------------------------------------------------------------

Ht0      =  37   ; % Ematocrito iniziale
W0       =  61   ; % [Kg] 
type_hdf = 'post'; % 'pre '  'post'
Vuf      =  1.3  ; % [Litri] Volume da ultrafiltrare
Td       =  240  ; % [min] tempo di dialisi
Qb       =  250  ; % [mL/min] portata ematica
Qs       =   90  ; % [mL/min] portata diluizione
Qd       =  375  ; % [mL/min] flusso dializzante

sacca_m = 'AX03' ; % codice sacca 'AX03' 'A161' sul monitor
sacca   = 'A161' ; % codice sacca 'AX03' 'A161' realmente usata
C_Na    = 140;     % [mmol/L] sodio prescritto
C_Bic   = 34;      % [mmol/L] HCO3- prescritto

eval(['save ',char(ID)])