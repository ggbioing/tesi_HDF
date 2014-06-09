function pre_processing(ID)

eval(['load ',char(ID)]); % carica i dati

n_sample = length(t_sample);
Cpl = zeros(n_sample,8);

soluto = {  'Na',  'K', 'Cl', 'Ca', 'P', 'Mg', 'urea', 'creat' };
MW     = [   23    39    35    40   96    24     60     113 ]; %peso molec.

%--------------------------------------------------------------------------
%-- TRAFORMAZIONE in mmol/L dei dati ematici per i soluti di interesse ----
%--------------------------------------------------------------------------
for i=1:n_sample   
    Cpl(i,:) = [ Na(i)  K(i)  Cl(i)  Ca(i)/MW(4)*10  P(i)/MW(5)*10 ...
                      Mg(i)/MW(6)*10  urea(i)/MW(7)*10  creat(i)/MW(8)*10 ];
end
Cpl0 = Cpl(1,:); % [mmol/L] Cpl dei soluti a inizio seduta
Tp0  = Tp(1);    % [gr/dL]  Proteine a inizio seduta

Ht0 = Ht0/100;    % [adim.]
Qb  = Qb/1000/60; % [L/sec]
Qs  = Qs/1000/60; % [L/sec]
Qd  = Qd/1000/60; % [L/sec]

%- Assegnare i parametri rho, k, eta se non sono ancora stati specificati -
if exist('RHO')==0, RHO=1; end
if exist('PAR')==0
    %           'Na' ,   'K' ,   'Cl',  'Ca' ,  'P' ,  'Mg' ,'urea', 'creat'
    PAR(:,1)=[ 2.5e-3 1.67e-4  1.67e-4 1.67e-4 3.5e-4 1.67e-4 1.3e-3  1.3e-4];
    PAR(:,2)=0.5*ones(8,1);
end %----------------------------------------------------------------------

eval(['save ',char(ID)])