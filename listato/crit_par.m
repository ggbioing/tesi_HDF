function J = J_CF(par)

global RHO PAR S paziente ...
    
eval(['load ',char(paziente)])

... % qui calcoli intermedi
    
if S==17
    rho = par;
else
    PAR(S)=par;
end

... % qui calcoli intermedi
    
%--------------------------------------------------------------------------
%- INTEGRAZIONE NUMERICA DEL SISTEMA DI ODE -------------------------------
%--------------------------------------------------------------------------
[TT,YY] = ode15s(@ODE_emodiafilt,[0:60:60*Td],Y0);

%--------------------------------------------------------------------------
%- RISULTATI --------------------------------------------------------------
%--------------------------------------------------------------------------
    MIC = YY(t_sample+1,[1:8])  ;  %  1..8  : Na..Urea intracellulare
    MEX = YY(t_sample+1,[1:8]+8);  % 12..19 : Na..Urea extracellulare
    
    VIC = YY(t_sample+1,17); % Volume intracellulare
    VIS = YY(t_sample+1,18); % Volume interstizio
    VPL = YY(t_sample+1,19); % Volume plasma
%_________________________________________________________________________|
TP = Tp0*YY(1,19)./VPL;
ad = Donnan(TP');

VTOT=VPL+VIS+VIC;
for i=1:length(t_sample)
    CIC(i,:) = MIC(i,:)./VIC(i);
    CIS(i,:) = MEX(i,:)./(VIS(i)+VPL(i)./ad([1:8],i)');
    SIM(i,:) = CIS(i,:)./ad([1:8],i)'; %CPL simulazione
end

eval(['load pesi_',char(paziente)])  % carica i pesi relativi a paz./seduta

E   = ([SIM TP] - [Cpl(:,[1:8]) Tp'])./[Cpl(:,[1:8]) Tp'];
E   = pesi(S,:)*E';
J   = sqrt( sum( E.^2 ) );
disp([num2str(S),' : ',num2str(par),' J = ',num2str(J)])
