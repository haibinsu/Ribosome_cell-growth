function [package, dist] = proteincombo(kcatKMc, kcatKMnc, qnc)

%to test different combination of rates at each Mg2+
%--> check effect of variation 

%concentration distribution
syms R Actc Actnc Act_starc Act_starnc 

Rtotal = 100;
T3c = 2;
T3nc = 15;

khyd = 500;
kpepc = 7;
kpepnc = 0.3;
qc = 1;  %assume constant , actually reducing ,as long as << kpepc 

Actcsol_store = NaN*ones(length(kcatKMc)*length(kcatKMnc)*length(qnc),1);
Actncsol_store = NaN*ones(length(kcatKMc)*length(kcatKMnc)*length(qnc),1);
Rsol_store = NaN*ones(length(kcatKMc)*length(kcatKMnc)*length(qnc),1);
Act_starcsol_store = NaN*ones(length(kcatKMc)*length(kcatKMnc)*length(qnc),1);
Act_starncsol_store = NaN*ones(length(kcatKMc)*length(kcatKMnc)*length(qnc),1);

count = 1;
for i = 1 : length(kcatKMc)
    for j = 1 : length(qnc)
        for k = 1 : length(kcatKMnc)
            %equation at steady state 
            eq1 = R * kcatKMc(i) * T3c - (khyd)*Actc == 0;
            eq2 = Actc * khyd - (qc + kpepc)*Act_starc == 0;

            eq3 =  R * kcatKMnc(k) * T3nc - (khyd)*Actnc == 0;
            eq4 = Actnc * khyd - (qnc(j) + kpepnc)*Act_starnc == 0;

            %eq5 =  (qc + kpepc)*Act_starc + (qnc + kpepnc)*Act_starnc - R * (k11 + k12) == 0

            eq7 = R + Actc + Actnc + Act_starc + Act_starnc == Rtotal ;

            eqn = [eq1, eq2, eq3, eq4, eq7];

            [Rsol, Actcsol, Actncsol, Act_starcsol, Act_starncsol ] = solve(eqn, [R Actc Actnc Act_starc Act_starnc]);
            Rsol_store(count,1) = Rsol;
            Actcsol_store(count,1) = Actcsol;
            Actncsol_store(count,1) = Actncsol;
            Act_starcsol_store(count,1) = Act_starcsol; 
            Act_starncsol_store(count,1) = Act_starncsol;
            
            Jtotal(count,1) = (Act_starcsol_store(count,1) *kpepc + Act_starncsol_store(count,1) *kpepnc)/Rtotal;
            Atotal(count,1) = Act_starcsol_store(count,1) *kpepc./(Act_starncsol_store(count,1) *kpepnc);
            JIS(count,1) =  (Actcsol_store(count,1)  + Actncsol_store(count,1) )*khyd/Rtotal;
            Q(count,1) = (Act_starcsol_store(count,1) *qc + Act_starncsol_store(count,1) .*qnc(j))/Rtotal;
            AIS(count,1) = Actcsol_store(count,1) ./Actncsol_store(count,1) ;
            Aproof(count,1) = (kpepc/(kpepc+qc))./(kpepnc./(kpepnc+qnc(j)));
            Jtotalc(count,1) = Act_starcsol_store(count,1) *kpepc/Rtotal;
            Jtotalnc(count,1) = Act_starncsol_store(count,1) *kpepnc/Rtotal;
            Qc(count,1) = Act_starcsol_store(count,1) *qc/Rtotal;
            Qnc(count,1) = Act_starncsol_store(count,1) .*qnc(j)/Rtotal;
            JISc(count,1) =  Actcsol_store(count,1) *khyd/Rtotal;
            JISnc(count,1) =  Actncsol_store(count,1) *khyd/Rtotal;    
            count = count + 1;
        end
    end
end

package = {Jtotal;Atotal;JIS;Q;AIS;Aproof;Jtotalc;Jtotalnc;Qc;Qnc;JISc;JISnc};
dist = {Rsol_store, Actcsol_store, Actncsol_store, Act_starcsol_store, Act_starncsol_store};
end