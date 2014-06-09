function dy = ODE_emodiafilt(t,y)

global Donnan hind rho Vrc Tp0 Tpis0 Vpl0 Vis0 Pac0 Pis0 alpha Rp gamma ...
       kf sigma k beta Gic type_hdf eff Td Quf Qf Qb  CL Cd Cs Qs Pvc   ...
       LP La Lv Pn_eq Cc Eis soluto
   
%--------------------------------------------------------------------------
%- VALORI ODE IN ENTRATA --------------------------------------------------
%--------------------------------------------------------------------------
for i=1:8                                                                 
    Mic(i) = y(i);   %  1..8  : Na..Gluc intracellulare
    Mex(i) = y(i+8); % 12..19 : Na..Gluc extracellulare
end    
    Vic = y(17); % Volume intracellulare
    Vis = y(18); % Volume plasmatico
    Vpl = y(19); % Volume interstiziale
    Pac = y(20); % Pressione ai capillari arteriosi
    Pis = y(21); % Pressione interstiziale
%--------------------------------------------------------------------------

Ht   = 1/(1+Vpl/Vrc);  % ematocrito
Tp   = Tp0*Vpl0/Vpl;   % proteinocrito
Tpis = Tpis0*Vis0/Vis; % proteine interstiziali
Fp   = Rp(Tp);         % frazione di acqua nel plasma
ad   = Donnan(Tp);     % coeff. di Donnan

%- Coefficiente di donnan ai capillari del dializzatore -------------------
if type_hdf == 'pre '
    ad_f = Donnan(Tp/(1+Qs/(Qb*(1-Ht))));
else
    ad_f=ad;
end
%--------------------------------------------------------------------------

%- Valori asintotici di equilibrio ( *_eq ) -------------------------------
Vtot_eq = Vic+Vis+Vpl ;
Vpl_eq  = 3/42*Vtot_eq;
ad_eq   = Donnan(Tp0*Vpl0/Vpl_eq);
Vis_eq  = 11/42*Vtot_eq ;
Vic_eq  = 28/42*Vtot_eq ;
Pac_eq  = Pac0 + 1/Cc*(Vpl_eq - Vpl0);
Pis_eq  = Pis0 + Eis*(Vis_eq-Vis0)   ;
Pn_eq   = ( (Pac_eq+Pvc)/2 - Pis_eq) ...
                         - ( LP(Tp0*Vpl0/Vpl_eq) - LP(Tpis0*Vis0/Vis_eq) );
Qic_eq = 0 ;
for i=1:8
    Mtot(i)   = Mic(i)+Mex(i);
    Cpl_eq(i) = Mtot(i)/(Vpl_eq + ad_eq(i)*(Vis_eq+beta(i)*Vic_eq));
    Cis_eq(i) = ad_eq(i)*Cpl_eq(i);
    Cic_eq(i) = beta(i)*Cis_eq(i) ;
    Qic_eq    = Qic_eq + gamma * kf * sigma(i) * ( Cic_eq(i) - Cis_eq(i) );
end
%--------------------------------------------------------------------------

Qic = 0 ;  % flusso intracellulare di liquido
for s=1:8  % Concentrazioni intra-cellulari, interstiz., plasm.  Na+...Gluc
    Cic(s) =  Mic(s)/Vic;
    Cis(s) =  Mex(s)/(Vpl/ad(s)+Vis);
    Cpl(s) =  Cis(s)/ad(s);
    Qic = Qic + gamma * kf * sigma(s) * ( Cic(s) - Cis(s) );
end
Qic_eff = (Qic - Qic_eq);

dVic = Qic_eff; % flusso di liquidi in entrata nel compart. intra-cellulare
for s=1:8   % Calcolo dei flussi di massa verso il compart. intra-cellulare
    phi(s)  = -  k(s)*( Cic(s) - beta(s)*Cis(s) );
    dMic(s) = phi(s) + Gic(s);
end

%- Calcolo della variazione della massa extracellulare --------------------
if type_hdf == 'pre '
    Qi = Qs + Qb*(1-Ht)*Fp;
    for s=1:8
       Ci(s)    = ( Qb*(1-Ht)*Cpl(s) + Qs*Cs(s) )/Qi;
       phi_t(s) = Qf*hind(s)*( Ci(s)+Cd(s) )/2;
       PHI(s)   = ( 1-Qf/Qi )*( CL(s)*( ad_f(s)*Ci(s)-Cd(s) )+phi_t(s) )...
                                          + Qf/Qi*(Qi-(1-eff(s))*Qf)*Ci(s);
       if PHI(s) > Qi*Ci(s), PHI(s) = Qi*Ci(s); end
       dMex(s)  = - phi(s) + Qs*Cs(s) - PHI(s);
       if t > Td*60, dMex(s)  = - phi(s); end
    end
elseif type_hdf == 'post'
    Qi = Qb*(1-Ht)*Fp;
    for s=1:8
       Ci(s)    = Cpl(s)/Fp;
       phi_t(s) = Qf*hind(s)*( Ci(s)+Cd(s) )/2;
       PHI(s)   = ( 1-Qf/Qi )*( CL(s)*( ad_f(s)*Ci(s)-Cd(s) )+phi_t(s) )...
                                          + Qf/Qi*(Qi-(1-eff(s))*Qf)*Ci(s);
       if PHI(s) > Qi*Ci(s), PHI(s)=Qi*Ci(s); end
       dMex(s)  = - phi(s) - PHI(s) + Qs*Cs(s);
       if t > Td*60, dMex(s)  = - phi(s); end
    end
end
%--------------------------------------------------------------------------

Pn = ( (Pac+Pvc)/2 - Pis )-( LP(Tp)-LP(Tpis) ); % Pres. di filtr. capillare

if t > Td*60, Quf=0; end    
if Vpl > 0 
    dVpl = - Quf - rho*(La+Lv)*(Pn - Pn_eq);
else % il volume plasmatico non può diventare negativo!
    dVpl = 0;
end

if Vis > 0 
    dVis = - Qic_eff + rho*(La+Lv)*(Pn - Pn_eq);
else % il volume interstiziale non può diventare negativo!
    dVis = 0;
end

if Pac>Pvc || dVpl > 0   
    dPac = 1/Cc * dVpl;    % variazione di pressione ai capillari arteriosi
else             % fenomeno del blocco di portata in un condotto non rigido
    dPac = 0;
end

dPis = Eis * dVis;     % variazione di pressione nel compart. interstiziale

%--------------------------------------------------------------------------
%- VALORI ODE IN USCITA ---------------------------------------------------
%--------------------------------------------------------------------------
for i=1:8
    dy(i)   = dMic(i);
    dy(i+8) = dMex(i);
end    
    dy(17) = dVic;
    dy(18) = dVis;
    dy(19) = dVpl;
    dy(20) = dPac;
    dy(21) = dPis;
%--------------------------------------------------------------------------
dy = dy';
