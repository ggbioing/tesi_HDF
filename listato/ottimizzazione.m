function ottimizzazione(ID)

eval([char(ID)])
pre_processing(char(ID))
global S RHO PAR paziente
paziente = ID;

for RIPETIZ=1:5    % ciclo per ripetere più volte l'ottimizz. dei parametri
    eval(['load ',char(paziente), ' RHO PAR'])
    P=[PAR(:); RHO];
    
    eval(['load ord_',char(ID),' y'])     % carica l'ordine di ottimizz.
    
    for S = y
        if S==17
            [P(S)] = fmincon( @J_CF,P(S),[],[],[],[], 0,10,[] );
        elseif S<=8
            [P(S)] = fmincon( @J_CF,P(S),[],[],[],[], 0, 1,[] );
        else 
            [P(S)] = fmincon( @J_CF,P(S),[],[],[],[],-1, 1,[] );
        end
        
        RHO    = P(end);
        PAR(:) = P(1:16);
        eval(['save ',char(paziente),' RHO PAR -append']) % aggiorna param.
    end
end