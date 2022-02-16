%include coupled nutrient and accuracy effect
%see how ksyneff - accuracy depend on [T3]
%assumption
kpepnc = 0.3; %1/s

khyd = 500; %1/s GTP hydrolysis is fast 
kpepc = 22; %1/s
qc = 1;  %cognate PR rejection rate - assume constant with Mg2+ ,as long as << kpepc 
%assumption taken from 10.1016/j.molcel.2005.12.018
Rtotal = 10; %uM
%total ternary complex = 100 uM 
T3c = 2; %uM 
T3nc = 15; %uM 

%for initial selection
%Johansson et al, PNAS, 2012 www.pnas.org/cgi/doi/10.1073/pnas.1116480109
kcatKMc_p1 = [60; 117; 147; 167; 180];  % cognate AAA uM^-1s^-1
kcatKMnc_p1 = [19; 66; 139; 327; 1750]; % near cognate GAA mM^-1s^-1 
kcatKMnc_p1 = kcatKMnc_p1/1000; %convert from mM^-1 to uM^-1

%for peptide bond
kcatKMpepnc_p1 = [3.9e-4; 2.7e-3; 9.86e-3;3.67e-2; 2.5e-1];  %unit is uM^-1s^-1
kcatKMpepc_p1 = [60; 117; 147; 167; 180]; %unit is uM^-1s^-1

%from coarse-grain scheme 
qnc_p1 = (kcatKMnc_p1./kcatKMpepnc_p1-1)*kpepnc; 


%try different set 
kcatKMc = [40; 60; 107; 135; 148; 170];  % cognate AAA uM^-1s^-1
kcatKMnc = [5; 19; 66; 119; 276; 1450]; % near cognate GAA mM^-1s^-1 
kcatKMnc = kcatKMnc/1000; %convert from mM^-1 to uM^-1

kcatKMpepnc = [1e-4; 3.9e-4; 2.4e-3; 4.86e-3;3.67e-2; 2.5e-1];  %unit is uM^-1s^-1
kcatKMpepc = [40; 60; 117; 147; 167; 180]; %unit is uM^-1s^-1

qnc = (kcatKMnc./kcatKMpepnc-1)*kpepnc; 


% kcatKMc = [40; 50; 60; 100; 140; 190];  % cognate AAA uM^-1s^-1
% qnc = [10; 8; 7; 4; 1.5; 0.8];
% kcatKMnc = [5;20; 40; 140; 500; 800]/1000; %convert from mM^-1 to uM^-1
% kcatKMpepnc = kcatKMnc./(qnc/kpepnc+1);
% 

%ksyneff expressed in Michaelis - Menten form 
fc = 0.02;
fnc = 0.15;
T3 = 100; %uM 
phiR = 0.2; 
NR = 12307; %aa/R

lambda = @(kel) phiR/NR*kel;
ksynmax = @(kfc,kpepc,kfnc,krejnc) (fc*kfc+fnc*kfnc./(1+krejnc./kpepnc))./(fc*kfc/kpepc+fnc*kfnc./krejnc);
KMeff = @(kfc,kpepc,kfnc,krejnc) 1./(fc*kfc/kpepc+fnc*kfnc./krejnc);
Accuracy = @(kfc,kpepc,kfnc,krejnc,krejc) fc*kfc./(fnc*kfnc).*(kpepc/(krejc+kpepc))./(kpepnc./(krejnc+kpepnc));
RPRc_numerator = @(kfc,kpepc,kfnc,krejnc,krejc) fc*kfc./(krejc+kpepc);  %remove T3 because it will be canceled in calculating cost
RPRnc_numerator = @(kfc,kpepc,kfnc,krejnc,krejc) fnc*kfnc./(krejnc+kpepnc);
speed_f = @(kfc,kpepc,kfnc,krejnc,T3) ksynmax(kfc,kpepc,kfnc,krejnc)*T3./(KMeff(kfc,kpepc,kfnc,krejnc)+T3);


cost_f = @(kfc,kpepc,kfnc,krejnc,krejc) (RPRc_numerator(kfc,kpepc,kfnc,krejnc,krejc)*krejc + RPRnc_numerator(kfc,kpepc,kfnc,krejnc,krejc).*krejnc)./(RPRc_numerator(kfc,kpepc,kfnc,krejnc,krejc)*kpepc + RPRnc_numerator(kfc,kpepc,kfnc,krejnc,krejc).*kpepnc);
cost_fnconly = @(kfc,kpepc,kfnc,krejnc,krejc) (RPRnc_numerator(kfc,kpepc,kfnc,krejnc,krejc).*krejnc )./(RPRc_numerator(kfc,kpepc,kfnc,krejnc,krejc)*kpepc + RPRnc_numerator(kfc,kpepc,kfnc,krejnc,krejc).*kpepnc);

ksynmaxlist_p1 = ksynmax(kcatKMc_p1,kpepc,kcatKMnc_p1,qnc_p1);
KMefflist_p1 = KMeff(kcatKMc_p1,kpepc,kcatKMnc_p1,qnc_p1);
Alist_p1 = Accuracy(kcatKMc_p1,kpepc,kcatKMnc_p1,qnc_p1,qc);

ksynmaxlist = ksynmax(kcatKMc,kpepc,kcatKMnc,qnc);
KMefflist = KMeff(kcatKMc,kpepc,kcatKMnc,qnc);
Alist = Accuracy(kcatKMc,kpepc,kcatKMnc,qnc,qc);
grlist = lambda(ksynmaxlist*T3./(T3+KMefflist)*3600);
costlist = cost_f(kcatKMc,kpepc,kcatKMnc,qnc,qc);



figure
yyaxis left
plot(Alist_p1, ksynmaxlist_p1, '.-')
hold on
scatter(Alist, ksynmaxlist)
yyaxis right
plot(Alist_p1, KMefflist_p1, '.-')
hold on
scatter(Alist, KMefflist)
set(gca,'XScale','log')

%physically possible range of ternary complex concentration
%high end 6e5 tRNA at 4 um3 --> 6e5/(6.02e23*4e-15)*1e6
%low end 2e4/(6.02e23*1e-15)*1e6 ~ 33 uM --> ~ 2000 ribosomes at 0.25 1/h
%with 2e4 tRNA and cell volume of 1 um3
T3list = [300 100 50 10];
figure
for i = 1 : length(T3list)
    T3 = T3list(i);
    plot(Alist, speed_f(kcatKMc,kpepc,kcatKMnc,qnc,T3),'.-','MarkerSize',10)
    hold on
end
legend('[T3] = 300 \muM','[T3] = 100 \muM', '[T3] = 50 \muM', '[T3] = 10 \muM')
set(gca,'XScale','log')
xlabel('Accuracy')
ylabel('k_{syn}^{eff} (1/s)')

%Accuracy recycling term
% a = 0.51;
% b = 2294;
% c = 0.48;
% % US mode
% T3_A = @(x,b,c) a*x./(b+x)+0.48;
c_prefactor = 1e-3;
b = 2000;
T3_A = @(x,b,c) 0.51./(1+exp(-c*(x-b)))+0.48;
figure
plot((1:1e4),T3_A(1:1e4,b,c_prefactor))
hold on
plot((1:1e4),T3_A(1:1e4,4000,c_prefactor))
plot((1:1e4),T3_A(1:1e4,6000,c_prefactor))
plot((1:1e4),T3_A(1:1e4,8000,c_prefactor))
plot((1:1e4),T3_A(1:1e4,1000,c_prefactor))

xlabel('Accuracy')
ylabel('Speed factor')

%different values of b depending on nutrient condition ?
% blist = [10000 5000 1800 800 100]; 
% blist = [5000 3000 1000 800 100]; 
%speed_recycling
speed_fA = @(kfc,kpepc,kfnc,krejnc,T3,A,b,c) speed_f(kfc,kpepc,kfnc,krejnc,T3).*T3_A(A,b,c);
figure
for i = 1 : length(T3list)
    T3 = T3list(i);
%     b = blist(i);
    plot(Alist, speed_fA(kcatKMc,kpepc,kcatKMnc,qnc,T3,Alist,b,c_prefactor),'.-','MarkerSize',10)
    hold on
end
legend('[T3] = 300 \muM','[T3] = 100 \muM', '[T3] = 50 \muM', '[T3] = 10 \muM')
set(gca,'XScale','log')
xlabel('Accuracy')
ylabel('k_{syn}^{eff} (1/s)')

figure
% plot(Alist,(fc*kcatKMc/kpepc)./(fnc*kcatKMnc./qnc),'sq-')
plot(Alist, ksynmaxlist,'k.-','MarkerSize',10)
hold on
scatter(Alist, fc*kcatKMc./(fc*kcatKMc/kpepc+fnc*kcatKMnc./qnc),'sq')
set(gca,'XScale','log')

figure
% plot(Alist,(fc*kcatKMc/kpepc)./(fnc*kcatKMnc./qnc),'sq-')
plot(Alist, KMefflist,'k.-','MarkerSize',20)
hold on
plot(Alist, 1./(fc*kcatKMc/kpepc),'ksq-')
set(gca,'XScale','log')

figure
plot(Alist, ksynmaxlist,'k.-','MarkerSize',20)
hold on

xlabel('Accuracy')
ylabel('k_{syn}^{max} (1/s)')
set(gca,'XScale','log')

figure
plot(Alist, KMefflist,'k.-','MarkerSize',20)
xlabel('Accuracy')
ylabel('K_{M}^{eff} (\muM)')
set(gca,'XScale','log')


figure
plot(Alist, ksynmaxlist./KMefflist,'k.-','MarkerSize',20)
xlabel('Accuracy')
ylabel('k_{syn}^{max}/K_{M}^{eff} (\muM^{-1}s^{-1})')
set(gca,'XScale','log')


figure
yyaxis left
plot(Alist, ksynmaxlist,'o-','MarkerFaceColor',[0, 0.4470, 0.7410])
xlabel('Accuracy')
ylabel('k_{syn}^{max} (1/s)')
ylim([10 24])
yyaxis right
plot(Alist, KMefflist,'sq-','MarkerFaceColor',[0.8500, 0.3250, 0.0980])
ylabel('K_{M}^{eff} (\muM)')
set(gca,'XScale','log')
xlim([10 1e5])

figure
plot(Alist, grlist,'k.-','MarkerSize',20)
xlabel('Accuracy')
ylabel('\lambda (1/h)')
set(gca,'XScale','log')

% figure
% % plot(Alist, costlist,'k.-','MarkerSize',20)
% % hold on
% bar(log10(Alist), [cost_fnconly(kcatKMc,kpepc,kcatKMnc,qnc,qc) costlist-cost_fnconly(kcatKMc,kpepc,kcatKMnc,qnc,qc)],'stack')
% legend('ncog','cog')
% xlabel('log(Accuracy)')
% ylabel('Cost')
% set(gca,'XScale','log')


% pick 5 points

% %analytical growth rate
% blist = [25000, 20000, 8000, 3000, 200];
% clist = [0.8e-3, 3e-3, 5e-3, 5e-3, 1e-2];
% m = @(x,b,c) 0.51./(1+exp(-c*(x-b)))+0.48;
% 
% aa_ana = @(kmet,ksynmax,KMeff,A,b,c) (NR*(KMeff*1e-6/c_factor)*kmet*Ke^2./(2*ksynmax.*m(A,b,c)*NP)).^(1/3);
% gr_ana = @(kmet,ksynmax,KMeff,A,b,c)  aa_ana(kmet,ksynmax,KMeff,A,b,c)./(3*NR*(KMeff*1e-6/c_factor)./(2*ksynmax.*m(A,b,c))+aa_ana(kmet,ksynmax,KMeff,A,b,c).*(NP/kmet+NR./(ksynmax.*m(A,b,c))));
% 
% lambda_ana = NaN*ones(length(kmetsample),length(Alist));
% 
% for i = 1 : length(kmetsample)
%     for j = 1 : length(Alist)
%         lambda_ana(i,j) = gr_ana(kmetsample(i),ksynmaxlist(i)*3600, KMefflist(i), Alist(j), blist(i), clist(i));
%     end
% end
% 
% figure
% for i = 1 : length(kmetsample)
%     plot(Alist, lambda_ana(i,:),'.-')
%     hold on
% end
% set(gca,'XScale','log')
% 
% figure
% for i = 1 : length(kmetsample)
%     [val, id] = max(lambda_ana(i,:));
%     scatter(Alist(id),aa_ana(kmetsample(i),ksynmaxlist(id)*3600, KMefflist(id), Alist(id), blist(i), clist(i)))
%     hold on
% end
%% relate competition effect on growth rate

%check competition_growthrate.m

%starting point is growth_HB_T3.m 

% remember to run bionumers.m

run('bionumers.m')

close all

% NR = 7336; %aa/ribosome
% NP = 300; %aa/metabolic protein 
% rho = 7.4e8; %aa/um^3 
% NA = 6.02e23;
% ksynmax = 22*3600; %aa/h
% lambdamax = 17.83*3600; %aa/h (synthesis speed in term of growth rate)
% % Ka = 2e-3; %M
% % Ke = 2e-1; %M
% RPpercent = 0.7;
% %dry weight = 1/3 wet weight
% dwconv = 0.33; 
% %new value of Ka Ke to try
% maa = 110*1.66e-24; % g = 110 Da %average mass of amino acid
%www.pnas.org/cgi/doi/10.1073/pnas.1421138111 

% ksynmax = 20*3600; %aa/h
Ka = 0.25e-3; %M
% Ke = 1e-2; %M %5e-2 1e-2
Ke = 2.5e-3;
Jmtm = 5.2e4/6 ; %convert ATP to aa because 6 ATP required for 1 peptide bond
%https://doi.org/10.1016/j.cels.2019.06.003

%use KMeff from part 1 for ternary complex instead of rescaling 
k2 = 100*3600; %1/h max charging speed 
% k2 = 80*3600;
R = 20; %uM 
% S = R/1.5; % --> charging synthase/ribosome ratio = 1/1.5 
S = R/(5); %S/R = 1/5 is more important 

c_factor = 0.0364; %constant factor between T3 and aa, get from growth_HB_T3.m 

%cell mass (aa)
M = @(y) NR*y(:,1) + NP*y(:,2) + y(:,3);

%volume (um^3)
V = @(y) M(y)/rho; %um^3 

%R mass fraction
Rmf = @(y) NR*y(:,1)./M(y);

%P mass fraction
Pmf = @(y) NP*y(:,2)./M(y);

%aa mass fraction
aamf = @(y) 1 - Rmf(y) - Pmf(y);


% %% use ksynmax and KMeff in part 2 - growth model
% Mpick = 0.8476e9; %aa
% nSim = 100;
% uRlist = 0.01:0.01:0.9;
% % Kaselect = 2.5e-4*KMefflist/KMefflist(3);
% Kaselect = 2.5e-4*ones(length(ksynmaxlist),1); 
% 
% % kmetsample = [0.07 0.1 0.14 0.17 0.24]'*3600;
% % kmetsample = [0.04 0.08 0.14 0.21 0.23]'*3600;
% % 
% % kmetsample = [0.01 0.05 0.18 0.25 0.3]'*3600;
% % nutrient = [1 3 5.3 7 11];
% % %assume kmetsample depends on nutrient as S shape
% % kmetmax = 0.3*3600;
% % kqual = 0.75;
% % kmet_N = @(n) kmetmax./(1+exp(-kqual*(n-5)));
% % n_list = (0:0.01:15)'; %nutrient quality
% % figure
% % plot(n_list,kmet_N(n_list)/3600)
% % hold on
% % scatter(nutrient, kmetsample/3600,'filled')
% % xlabel('Nutrient quality')
% % ylabel('k_{met}^{eff} (1/s)')
% 
% kmetsample = [0.05 0.1 0.18 0.25 0.3]'*3600;
% nutrient = [0.7 2 5.3 7 11];
% %assume kmetsample depends on nutrient as S shape
% kmetmax = 0.3*3600;
% kqual = 0.5;
% kmet_N = @(n) kmetmax./(1+exp(-kqual*(n-4)));
% n_list = (0:0.01:15)'; %nutrient quality
% figure
% plot(n_list,kmet_N(n_list)/3600)
% hold on
% scatter(nutrient, kmetsample/3600,'filled')
% xlabel('Nutrient quality')
% ylabel('k_{met}^{eff} (1/s)')
% 
% Mx = 1.65e9;
% My = 2e8;
% M0f = @(x) My + Mx*x;
% % M0f(grmax_nutrient) 
% 
% Mpicklist = 1e9*[1.0091,    1.6177,    2.5050,    3.0691,    3.4578];
% % Mlist = 1e9*[0.976 1.807 2.531 2.791 2.951];    
%     
% %effect of accuracy on translation 
% 
% par2.NR = NR;
% par2.NP = NP;
% par2.rho = rho;
% par2.c_factor = c_factor;
% par2.NA = NA;
% par2.Mpick = Mpick;
% par2.nSim = nSim;
% par2.Ke = Ke;
% par2.b = b;
% par2.c = c_prefactor;
% 
% odestore_kmet = cell(length(kmetsample),1);
% odemax_kmet = cell(length(kmetsample),1);
% kelongstore_kmet = cell(length(kmetsample),1);
% Alist_kmet = cell(length(kmetsample),1);
% lambda_ana_kmet = cell(length(kmetsample),1); 
% koriginal_kmet = cell(length(kmetsample),1);
% flux_kmet = cell(length(kmetsample),1);
% 
% % blist = [1000, 2000, 3000, 4000, 5000];
% % blist = [6000, 4000, 1000, 500, 200];
% 
% %use this
% blist = [7000, 6000, 2000, 1000, 500];
% clist = [1e-3, 1e-3, 1e-3, 1e-3, 5e-2];
% T3_A = @(x,b,c) 0.51./(1+exp(-c*(x-b)))+0.48;
% 
% color_list = [77 191 237; 113 152 241; 148 114 244; 188 71 248; 255 0 255]/255; %list of color
% 
% 
% 
% 
% % degrade = @(x) 1*exp(-1e-3*(x-1000))./(10+exp(-5e-4*(x-1000)));
% % degrade = @(x) 0.1*exp(-1e-3*(x-1000))./(1-0.1+0.1*exp(-1e-3*(x-1000)));
% degrade = @(x) 0.3./(1 + 2*exp(x./1000));
% % degrade = @(x) 20./x;
% plot((10:1e5), degrade((10:1e5)))
% set(gca,'XScale','log')
% xlabel('Accuracy')
% ylabel('degradation rate (1/h)')
% 
% figure
% for i = 1 : length(kmetsample)
%     plot((10:1e5), T3_A((10:1e5),blist(i),clist(i)), 'Color', color_list(i,:))
%     hold on    
% end
% set(gca,'XScale','log')
% xlabel('Accuracy')
% ylabel('m(A)')
% 
% 
% for m = 1: length(kmetsample)
%     kmet = kmetsample(m);
%     par2.kmet = kmet;
%     par2.b = blist(m);
%     par2.c = clist(m);
%     par2.Mpick = Mpicklist(m);
% %     poolobj = gcp;
% %     addAttachedFiles(poolobj,{'ode_growth_linking.m'})
%     %
%     [odestore_kmet{m}, odemax_kmet{m}, kelongstore_kmet{m}, koriginal_kmet{m}, Alist_kmet{m}, lambda_ana_kmet{m}, flux_kmet{m}] = lambda_A_T3A(ksynmaxlist,KMefflist,Kaselect,Alist,uRlist,M,par2,T3_A,degrade);
%     
% end
%   
% 
% 
% % %analytical growth rate
% % qA = @(x,b,c) 1;
% % 
% % c_prefactor = 1e-3;
% % % c_prefactor = 100;
% % b = 2000;
% % qA = @(x,b,c) 0.51./(1+exp(-c*(x-b)))+0.48;
% % blist = [7000, 6000, 2000, 1000, 500];
% % clist = [1e-3, 1e-3, 1e-3, 1e-3, 5e-2];
% % for i = 1 : length(blist)    
% %     plot((1:1e4), qA((1:1e4),blist(i),clist(i)))
% %     hold on
% % end
% % aa_ana = @(ksynmax,KMeff,A) (NR*(KMeff*1e-6/c_factor)*kmet*Ke^2./(2*ksynmax.*qA(A)*NP)).^(1/3);
% % gr_ana = @(ksynmax,KMeff,A)  aa_ana(ksynmax,KMeff,A)./(3*NR*(KMeff*1e-6/c_factor)./(2*ksynmax.*qA(A))+aa_ana(ksynmax,KMeff,A).*(NP/kmet+NR./(ksynmax.*qA(A))));
% % 
% % for i = 1: length(kmetsample)
% %     kmet = kmetsample(i);
% %     b = blist(i);
% %     c_prefactor = clist(i);
% %     aa_ana = @(ksynmax,KMeff,A,b,c) (NR*(KMeff*1e-6/c_factor)*kmet*Ke^2./(2*ksynmax.*qA(A,b,c)*NP)).^(1/3);
% %     gr_ana = @(ksynmax,KMeff,A,b,c)  aa_ana(ksynmax,KMeff,A,b,c)./(3*NR*(KMeff*1e-6/c_factor)./(2*ksynmax.*qA(A,b,c))+aa_ana(ksynmax,KMeff,A,b,c).*(NP/kmet+NR./(ksynmax.*qA(A,b,c))));
% %     lambda_ana = gr_ana(ksynmaxlist*3600, KMefflist, Alist, b, c_prefactor);
% %     plot(Alist, lambda_ana,'.-')
% %     hold on
% % end
% % set(gca,'XScale','log')
% % ylim([0 2])
% 
% 
% % c = [77 191 237; 113 152 241; 148 114 244; 188 71 248; 255 0 255]/255; %list of color
% 
% %mass
% % Mf = @(y) NR*y(1)+NP*y(2)+y(3);
% %volume
% % Vf = @(y) Mf(y)/rho;
% %elongation rate
% %c_factor to convert aa to Ternary complex
% kelongf = @(kmax,KM,y,Alist,b,c) T3_A(Alist,b,c).*3600.*kmax.*c_factor.*y(:,3)./(c_factor*y(:,3)+KM.*1e-6.*NA.*V(y)*1e-15);
% %metabolic rate
% kmetf = @(y,kmet) kmet./(1+(y(:,3)./(Ke*NA*V(y)*1e-15)).^2);
% %input flux
% JM_f = @(y,kmet) y(:,2).*kmetf(y,kmet);
% %output flux
% JS_f = @(y,kmax,KM,Alist,b,c) y(:,1).*kelongf(kmax,KM,y,Alist,b,c);
% 
% fluxstore = NaN*ones(length(Alist),length(kmetsample)*2);
% %1st -2nd col: input flux, growth rate flux for 1st kmet
% %3rd - 4th col: input flux, growth rate flux for 2nd kmet
% 
% JSstore = NaN*ones(length(Alist),length(kmetsample));
% %1st col: synthesis flux for 1st kmet, 
% %2nd col: synthesis flux for 2nd kmet,
% 
% pick_Q = NaN*ones(length(kmetsample),12);
% %1st col: intake flux, 2nd - 12th col: similar to odemax_kmet 
% 
% pick_Sd = NaN*ones(length(kmetsample),3);
% %1st col: synthesis flux, 2nd col: degradation flux, 3rd col: lambda*M
% 
% for i = 1 : length(kmetsample)
%    fluxstore(:,2*i-1:2*i) = [JM_f(odemax_kmet{i},kmetsample(i)),  odemax_kmet{i}(:,4)*Mpicklist(i)];
%    JSstore(:,i) = JS_f(odemax_kmet{i},ksynmaxlist,KMefflist,Alist,blist(i),clist(i));
%    [val, id] = max(odemax_kmet{i}(:,4));
%    pick_Q(i,:) = [fluxstore(id,2*i-1), odemax_kmet{i}(id,:)];
%    pick_Sd(i,:) = [JSstore(id,i), flux_kmet{m}(id,3), odemax_kmet{i}(id,4)*Mpicklist(i)];
% end
% 
% %check flux balance for each nutrient condition (i.e kmet)
% m = 1;
% %check elongation rate 
% kelongf(ksynmaxlist,KMefflist,odemax_kmet{m},Alist,blist(m),clist(m))./kelongstore_kmet{m}
% 
% 
% %check flux balance
% %input flux check: numerical vs analytical 
% flux_kmet{m}(:,1)./fluxstore(:,2*m-1)
% %synthesis flux check
% flux_kmet{m}(:,2)./JSstore(:,m)
% %degradation flux check
% flux_kmet{m}(:,3)./(degrade(Alist).*(NR*odemax_kmet{m}(:,1)+NP*odemax_kmet{m}(:,2)))
% % R*kelong + lambda*aa = JM
% (flux_kmet{m}(:,2)+odemax_kmet{m}(:,4).*odemax_kmet{m}(:,3))./flux_kmet{m}(:,1)
% % JM - degrade flux = lambda * M
% (flux_kmet{m}(:,3)+odemax_kmet{m}(:,4).*Mpicklist(m))./flux_kmet{m}(:,1)
% 
% 
% namelist = {'R','P','aa'};
% masslist = [NR NP 1];
% 
% for i = 1 : length(namelist)
%     
%     yL = sprintf('%s^* (molecules/cell)',namelist{i}); 
%     yR = ['\phi_' sprintf('{%s}^*',namelist{i})];       
%     figure
%     yyaxis left
%     plot( pick_Q(:,5),   pick_Q(:,i+1),'.-','MarkerSize',20)
%     ylabel(yL)
%     xlabel('Growth rate (1/h)')
%     yyaxis right
%     plot(pick_Q(:,5),  pick_Q(:,i+7),'sq-')
%     ylabel(yR)        
% end
% 
% % figure
% % plot(pick_Q(:,5), pick_Q(:,1), 'sq-', 'MarkerSize',10)
% % hold on
% % scatter(pick_Q(:,5), pick_Q(:,5).*Mpicklist')
% 
% figure
% yyaxis left
% bar(pick_Q(:,5), [pick_Sd(:,3), pick_Sd(:,2)], 'stacked')
% hold on
% scatter(pick_Q(:,5), pick_Q(:,1), 'filled','m')
% xlabel('Growth rate (1/h)')
% ylabel('Flux (aa/h)')
% % legend('Numerical','Analytical')
% yyaxis right
% plot(pick_Q(:,5), Agrmax_store, 'sq-','MarkerSize',7,'MarkerFaceColor', [0.85 0.33 0.1])
% set(gca,'YScale','log')
% ylabel('Accuracy')
% 
% legend('Growth flux','Degradation flux','Input flux','Accuracy')
% 
% 
% grmax_store = NaN*ones(length(kmetsample),1);
% Agrmax_store = NaN*ones(length(kmetsample),1);
% figure
% for m = 1 : length(kmetsample)    
%     plot(Alist, odemax_kmet{m}(:,4),'.-','Color',color_list(m,:))
% %     plot(nutrient, odemax_kmet{m}(:,4),'.-','Color',c(m,:))
% 
%     hold on
%     [val, id] = max(odemax_kmet{m}(:,4));
%     grmax_store(m) = val;
%     Agrmax_store(m) = Alist(id);
%     scatter(Agrmax_store(m), grmax_store(m), [],color_list(m,:),'filled')
% %     scatter(nutrient(id), val, [],c(m,:),'filled')
% 
%     %analytical solution 
% %     scatter(Alist, lambda_ana_kmet{m},[],c(m,:),'filled')
% end
% set(gca,'XScale','log')
% xlabel('Accuracy')
% % xlabel('Nutrient quality')
% ylabel('Growth rate (1/h)')
% 
% figure
% plot(pick_Q(:,5),pick_Q(:,5)./pick_Q(:,1), 'sq-', 'MarkerSize',5)
% xlabel('Growth rate (1/h)')
% ylabel('Growth rate/Input flux (1/aa)')
% 
% 
% 
% figure
% plot(Agrmax_store,pick_Q(:,5)./pick_Q(:,1), 'sq-', 'MarkerSize',5)
% xlabel('Accuracy')
% ylabel('Growth rate/Input flux (1/aa)')
% set(gca,'XScale','log')
% 
% 
% figure
% for i = 1 : length(kmetsample)
% %    plot(Alist,  fluxstore(:,2*i)./fluxstore(:,2*i-1),'.-')
%    plot(Alist,  fluxstore(:,2*i),'.-','Color',color_list(i,:))
%    hold on
%    scatter(Agrmax_store(i),pick_Q(i,1),[],color_list(i,:),'filled')
% end
% set(gca,'XScale','log')
% xlabel('Accuracy')
% ylabel('Input flux (aa/h)')
% 
% figure
% yyaxis left 
% plot(nutrient, grmax_store, '.-')
% hold on
% scatter(nutrient, grmax_store,'filled')
% ylabel('Growth rate (1/h)')
% yyaxis right
% plot(nutrient, Agrmax_store, '.-')
% hold on
% scatter(nutrient, Agrmax_store,'filled')
% xlabel('Nutrient quality')
% ylabel('Accuracy')
% set(gca,'YScale','log')
% 
% figure
% for m = 1 : length(kmetsample)
%     plot(Alist, kelongstore_kmet{m}/3600,'.-','Color',color_list(m,:))
%     hold on
% end 
% set(gca,'XScale','log')
% xlabel('Accuracy')
% ylabel('k^{eff}_{syn} (1/s)')
% 
% figure
% for m = 1 : length(kmetsample)
%     plot(Alist, koriginal_kmet{m}/3600,'.-','Color',color_list(m,:))
%     hold on
% end 
% set(gca,'XScale','log')
% xlabel('Accuracy')
% ylabel('k_{el} (1/s)')
% 
% %hand picked for China
% idx_A = [5,4,3,2];
% China_data = NaN*ones(length(idx_A),11);
% China_flux = NaN*ones(length(idx_A),3);
% %1st - 2nd - 3rd col: Accuracy - input flux - flux for growth (lambda*M)
% figure
% for j = 1 : length(kmetsample)
%    plot(Alist, odemax_kmet{j}(:,4),'.-','Color',color_list(j,:)) 
%    hold on 
%    if j < 5
%    scatter(Alist(idx_A(j)), odemax_kmet{j}(idx_A(j),4),[],color_list(j,:),'filled')
%    China_data(j,:) = odemax_kmet{j}(idx_A(j),:);
%    China_flux(j,:) = [Alist(idx_A(j)), fluxstore(idx_A(j),2*j-1:2*j)];
%    end
% end
% set(gca,'XScale','log')
% xlabel('Accuracy')
% ylabel('Growth rate (1/h)')
% 
% figure
% plot(China_flux(:,1), China_flux(:,2),'.-','MarkerSize',20)
% hold on
% plot(China_flux(:,1), China_flux(:,3),'.-','MarkerSize',20)
% set(gca,'XScale','log','YScale','log')
% xlabel('Accuracy')
% ylabel('Flux/cell (aa/h)')
% legend('Input flux','Flux for growth')
% 
% %hand picked for US
% idx_A = [1, 2, 3, 3];
% US_data = NaN*ones(length(idx_A),11);
% US_flux = NaN*ones(length(idx_A),3);
% %1st - 2nd - 3rd col: Accuracy - input flux - flux for growth (lambda*M)
% 
% figure
% for j = 1 : 5
%    plot(Alist, odemax_kmet{j}(:,4),'.-','Color',color_list(j,:)) 
%    hold on 
%    if j < 5
%    scatter(Alist(idx_A(j)), odemax_kmet{j}(idx_A(j),4),[],color_list(j,:),'filled')
%    US_data(j,:) = odemax_kmet{j}(idx_A(j),:);
%    US_flux(j,:) = [Alist(idx_A(j)), fluxstore(idx_A(j),2*j-1:2*j)];
% 
%    end
% end
% set(gca,'XScale','log')
% xlabel('Accuracy')
% ylabel('Growth rate (1/h)')
% 
% figure
% plot(US_flux(:,1), US_flux(:,2),'.-','MarkerSize',20)
% hold on
% plot(US_flux(:,1), US_flux(:,3),'.-','MarkerSize',20)
% set(gca,'XScale','log')
% xlabel('Accuracy')
% ylabel('Flux/cell (aa/h)')
% legend('Input flux','Flux for growth')
% 
% figure
% plot(US_flux(:,2),US_data(:,4),'.-','MarkerSize',20)
% hold on
% plot(China_flux(:,2),China_data(:,4),'.-','MarkerSize',20)
% ylabel('Growth rate (1/h)')
% xlabel('Input flux/cell (aa/h)')
% legend('US','China')
% 
% %plot Ribosome
% namelist = {'R','P','aa'};
% masslist = [NR NP 1];
% 
% mode_gr = 'China';
% mode_gr = 'US';
% data_plot = NaN*ones(length(idx_A),11);
% if strcmp(mode_gr,'China')
%     %China
%     idx_A = [5,4,3,2];
%     data_plot = China_data;
% else
%     %US
%     idx_A = [1, 2, 3, 3];
%     data_plot = US_data;
% end
% for i = 1 : length(namelist)
%     
%     yL = sprintf('%s^* (molecules/cell)',namelist{i}); 
%     yR = ['\phi_' sprintf('{%s}^*',namelist{i})];       
%     figure
%     yyaxis left
%     plot( data_plot(:,4),  data_plot(:,i),'.-','MarkerSize',20)
%     ylabel(yL)
%     xlabel('Growth rate (1/h)')
%     yyaxis right
%     plot(data_plot(:,4), data_plot(:,i+6),'sq-')
%     ylabel(yR)    
%     title(mode_gr)
% end
% 

%% explicit simulation instead of putting in function

Mpick = 0.8476e9; %aa
nSim = 100;
uRlist = 0.01:0.01:0.9;
% Kaselect = 2.5e-4*KMefflist/KMefflist(3);
Kaselect = 2.5e-4*ones(length(ksynmaxlist),1); 

% kmetsample = [0.07 0.1 0.14 0.17 0.24]'*3600;
% kmetsample = [0.04 0.08 0.14 0.21 0.23]'*3600;
% 
% kmetsample = [0.01 0.05 0.18 0.25 0.3]'*3600;
% nutrient = [1 3 5.3 7 11];
% %assume kmetsample depends on nutrient as S shape
% kmetmax = 0.3*3600;
% kqual = 0.75;
% kmet_N = @(n) kmetmax./(1+exp(-kqual*(n-5)));
% n_list = (0:0.01:15)'; %nutrient quality
% figure
% plot(n_list,kmet_N(n_list)/3600)
% hold on
% scatter(nutrient, kmetsample/3600,'filled')
% xlabel('Nutrient quality')
% ylabel('k_{met}^{eff} (1/s)')

%use this
kmetsample = [0.05 0.1 0.18 0.25 0.3]'*3600;
nutrient = [0.7 2 5.3 7 11];

% kmetsample = [0.05 0.1 0.14 0.18 0.22]'*3600;
% nutrient = [0.7 2 4.3 6.4 11];

%assume kmetsample depends on nutrient as S shape
kmetmax = 0.3*3600;
kqual = 0.5;
kmet_N = @(n) kmetmax./(1+exp(-kqual*(n-4)));
n_list = (0:0.01:15)'; %nutrient quality
figure
plot(n_list,kmet_N(n_list)/3600)
hold on
scatter(nutrient, kmetsample/3600,'filled')
xlabel('Nutrient quality')
ylabel('k_{met}^{eff} (1/s)')

Mx = 1.65e9;
My = 2e8;
M0f = @(x) My + Mx*x;
% M0f(grmax_nutrient) 

Mpicklist = 1e9*[0.9744,    1.5589,    2.3683,    3.0574,    3.4404];
% Mlist = 1e9*[0.976 1.807 2.531 2.791 2.951];    
    
% blist = [1000, 2000, 3000, 4000, 5000];
% blist = [6000, 4000, 1000, 500, 200];
blist = [7000, 6000, 2000, 1000, 500];
% clist = [1e-3, 1e-3, 1e-3, 1e-3, 5e-2];
clist = [1e-3, 1e-3, 1e-3, 1e-3, 1e-2];

color_list = [77 191 237; 113 152 241; 148 114 244; 188 71 248; 255 0 255]/255; %list of color
T3_A = @(x,b,c) 0.51./(1+exp(-c*(x-b)))+0.48;


%another blist due to extra point
% blist = [30000, 20000, 5000, 1000, 500];
% clist = [0.8e-3, 1e-3, 1e-3, 1e-3, 1e-2];

%not too bad
% blist = [25000, 20000, 8000, 5000, 500];
% clist = [0.8e-3, 3e-3, 5e-4, 5e-3, 1e-2];

%not too bad
% blist = [25000, 20000, 8000, 5000, 1000];
% clist = [0.8e-3, 3e-3, 5e-4, 2e-3, 5e-3];
% degrade = @(x) 0.3./(1 + 2*exp(x./1000));

%from 14 of test_package
blist = [25000, 18000, 5000, 1000, 500];
clist = [1e-4, 0.25e-3, 0.5e-3, 8e-4, 2e-3];
T3_A = @(x,b,c) 1.2./(2+exp(-c*(x-b)))+0.4;

degrade = @(x) 0.3./(1 + 2*exp(x./1000));

figure
plot((10:1e5), degrade((10:1e5)))
set(gca,'XScale','log')
xlabel('Accuracy')
ylabel('Degradation rate (1/h)')

figure
for i = 1 : length(kmetsample)
    plot((10:1e5), T3_A((10:1e5),blist(i),clist(i)), 'Color', color_list(i,:))
    hold on    
end
set(gca,'XScale','log')
xlabel('Accuracy')
ylabel('m(A)')

%effect of accuracy on translation 
BA = 1000; 
% qA = @(x) x./(x+BA);
% qA = @(x) 1 - exp(-x/BA);
% qA = @(x) exp(-0.008*NR./x);
qA = @(x) exp(-0.01*NR./x);  %x is accuracy 
    
odestore_kmet = cell(length(kmetsample),1);
odemax_kmet = cell(length(kmetsample),1);
kelongstore_kmet = cell(length(kmetsample),1);
korginalstore_kmet = cell(length(kmetsample),1);
Alist_kmet = cell(length(kmetsample),1);
lambda_ana_kmet = cell(length(kmetsample),1); 
koriginal_kmet = cell(length(kmetsample),1);
flux_kmet = cell(length(kmetsample),1);

%store the fastest growth rate and corresponding accuracy in each nutrient condition 
grmax_nutrient = NaN*ones(length(kmetsample),1);
ksynmax_nutrient = NaN*ones(length(kmetsample),1);
KM_nutrient = NaN*ones(length(kmetsample),1);
kelong_nutrient = NaN*ones(length(kmetsample),1);
koriginal_nutrient =  NaN*ones(length(kmetsample),1);
A_nutrient = NaN*ones(length(kmetsample),1);
flux_nutrient = NaN*ones(length(kmetsample),3);
%1st - 2nd - 3rd col: input flux - synthesis flux - degradation flux 
cost_nutrient = NaN*ones(length(kmetsample),1);
chosen = NaN*ones(length(kmetsample),11);

for m = 1 : length(kmetsample)
    kmet = kmetsample(m);
    Mth = Mpicklist(m);
    b = blist(m);
    c = clist(m);
    
    odestore = cell(length(ksynmaxlist),1);
    odemax = NaN*ones(length(ksynmaxlist),11);
    kelongstore = NaN*ones(length(ksynmaxlist),1);
    koriginal = NaN*ones(length(ksynmaxlist),1);
    fluxstore = NaN*ones(length(ksynmaxlist),1);
    %1st col: input flux, 2nd col: synthesis flux, 3rd col: degrade flux
    %each row corresponds to each ksynmax-KMeff-Accuracy
   
    for j = 1 : length(ksynmaxlist)
    
        ksynmaxT3 = ksynmaxlist(j)*3600; %1/h
        KM = KMefflist(j)*1e-6; %M
        KMeff = KM;
        Ka = Kaselect(j); %M 
        A = Alist(j);
        %cell mass (aa)
        Mf = @(y) NR*y(1)+NP*y(2)+y(3);
        %volume (um^3)
        Vf = @(y) Mf(y)/rho; %um^3 
        %R mass fraction
        Rmfc = @(y) NR*y(1)./Mf(y);

        %P mass fraction
        Pmfc = @(y) NP*y(2)./Mf(y);

        %aa mass fraction
        aamfc = @(y) 1 - Rmfc(y) - Pmfc(y);


    %     T3convertf_T3 = @(y) k2*S/R*KM*NA*Vf(y)*1e-15*y(3)./(ksynmax*Ka*NA*Vf(y)*1e-15+(ksynmax-k2*S/R)*y(3));
        %linear mapping aa to T3   
    %     T3convertf_T3 = @(y) k2*S/R*KM*y(3)./(ksynmaxT3*Ka);
        T3convertf_T3 = @(y) c_factor*y(3);
        %translation rate based on T3
%         ktransf_T3 = @(y) ksynmaxT3*T3convertf_T3(y)./(T3convertf_T3(y)+KMeff*NA*Vf(y)*1e-15);
    %include accuracy damping effect
%     if growth_mode == 1
        ktransf_T3 = @(y) ksynmaxT3*T3convertf_T3(y)./(T3convertf_T3(y)+KMeff*NA*Vf(y)*1e-15);
%     else
%         ktransf_T3 = @(y,A,b,c) ksynmaxT3*T3convertf_T3(y)*T3_A(A,b,c)./(T3convertf_T3(y)+KMeff*NA*Vf(y)*1e-15);
%     end
        %full mapping [aa] to [T3]
        %T3convertf_T3 gives [T3] concentration 
    %     T3convertf_T3 = @(y) k2*S/R*KM*y(3)./(ksynmaxT3*(Ka*NA*Vf(y)*1e-15+y(3))-S/R*k2*y(3));
    %     ktransf_T3 = @(y) ksynmaxT3*T3convertf_T3(y)./(T3convertf_T3(y)+KMeff);


    %     ktransf_T3 = @(y) ksynmaxT3*y(3)./(y(3)+Ka*NA*Vf(y)*1e-15);

        flist = {Mf,Vf,ktransf_T3,T3convertf_T3};

        %metabolic rate
        kmetf = @(y) kmet./(1+(y(3)/(Ke*NA*Vf(y)*1e-15))^2);
        flist{5} = kmetf; 
        flist{6} = qA; 
        flist{7} = T3_A;
        flist{8} = degrade;
        
        store = NaN*ones(length(uRlist),10);
%         Mth = Mpick;
        
        parfor i =  1 : length(uRlist)

            R0 = 2e4;
            P0 = 1e5;
            M0 = 3e8; %aa
            aa0 = M0 - (NR*R0 + P0*NP);
    %         R0 = nstore(j,1);
    %         P0 = nstore(j,2);
    %         aa0 = nstore(j,3);
            y0 = [R0 P0 aa0];
%             if growth_mode == 1  %US mode
%                 c = 1; %dummy
%                 par = [ksynmaxT3 kmet uRlist(i) NR NP A b c];
%             else
                par = [ksynmaxT3 kmet uRlist(i) NR NP A b c];
%             end
            k = 1;
            t_tot = 0;
            y_tot = y0;
            tspan = [0 200];  
            gr = 1;
            Opt = odeset('Events', @(t,y) myEvent_growth_opt2(t,y,Mth,par));
            %run to reach exponential steady state 
            while k <= nSim
                [t,y, te, ye, ie] = ode15s(@(t,y)  ode_growth_linkingT3A(t,y,par,flist),tspan,y0, Opt);
                if isempty(ye) ~= 1
                    t_tot = [t_tot ; [t;te] + t_tot(end)];
                    y_tot = [ y_tot ; [y;ye] ];    
                    y0 = ye/2;
                    k = k + 1;
    %                 dy = ode_growth_findmax(t(end),y(end,:),par);
    %                 gr = sum(dy.*[NR;NP;1])/sum(y(end,:).*[NR NP 1]);
                    %growth rate at current division sets future threshold mass 
    %                 Opt = odeset('Events', @(t,y) myEvent_growth_opt2(t,y,Mth,par));            
                else
                    gr = NaN;
                    break;
                end
            end
        %     semilogy(t_tot,y_tot)
            if isnan(gr) ~= 1 
            [t,y,te,ye,ie] = ode15s(@(t,y)  ode_growth_linkingT3A(t,y,par,flist),[0:1e-3:2*te],y0,Opt);
        %     semilogy(t,y)
            hfit = dy_dt_eval(t,par,flist,y,M);
            store(i,:) = [y0 hfit.p1 uRlist(i) 1 - uRlist(i) Rmfc(y0) Pmfc(y0) aamfc(y0) T3convertf_T3(y0)];

%             store(i,:) = [y0 log(2)/t(end) uRlist(i) 1 - uRlist(i) Rmfc(y0) Pmfc(y0) aamfc(y0) T3convertf_T3(y0)];
            

            end
            
          
        end
    %     plot(uRlist,store(:,4))
        odestore{j} = store;
        %pick uR parition that gives fastest growth rate 
        [val, idx] = max(store(:,4));
        odemax(j,:) = [store(idx,:) M(store(idx,:))];

        kelongstore(j) = ktransf_T3(odemax(j,:))*T3_A(A,b,c);
        koriginal(j) = ktransf_T3(odemax(j,:));
        
        fluxstore(j,1) = odemax(j,2)*kmetf(odemax(j,:));
        fluxstore(j,2) = odemax(j,1)*kelongstore(j);
        fluxstore(j,3) = degrade(A)*(NR*odemax(j,1)+NP*odemax(j,2));
    end 
%     for j = 1 : length(ksynmaxlist)
%         figure 
%         plot(odestore{j}(:,5),odestore{j}(:,4))
%     end
    %analytical growth rate
    aa_ana = @(ksynmax,KMeff,A,b,c) (NR*(KMeff*1e-6/c_factor)*kmet*Ke^2./(2*ksynmax.*qA(A).*T3_A(A,b,c)*NP)).^(1/3);
    gr_ana = @(ksynmax,KMeff,A,b,c)  aa_ana(ksynmax,KMeff,A,b,c)./(3*NR*(KMeff*1e-6/c_factor)./(2*ksynmax.*qA(A).*T3_A(A,b,c))+aa_ana(ksynmax,KMeff,A,b,c).*(NP/kmet+NR./(ksynmax.*qA(A).*T3_A(A,b,c))));

    lambda_ana = gr_ana(ksynmaxlist*3600, KMefflist, Alist, b, c);

%     scatter(Alist, lambda_ana)
%     hold on
%     plot(Alist, odemax(:,4))
%     
    odestore_kmet{m} = odestore; 
    odemax_kmet{m} = odemax; 
    kelongstore_kmet{m} = kelongstore;
    korginalstore_kmet{m} = koriginal;
    flux_kmet{m} = fluxstore;
    Alist_kmet{m} = Alist;
    lambda_ana_kmet{m} = lambda_ana;
    
    %store fastest growth rate and corresponding quantities for each nutrient condition
    %due to variation in accuracy, ksynmax and KMeff
    [val, id] = max(odemax(:,4));  %max growth rate
    grmax_nutrient(m) = val;
    ksynmax_nutrient(m) = ksynmaxlist(id);
    KM_nutrient(m) = KMefflist(id);
    kelong_nutrient(m) = kelongstore(id);
    koriginal_nutrient(m) = koriginal(id);
    A_nutrient(m) = Alist(id);
    flux_nutrient(m,:) = fluxstore(id,:);
    cost_nutrient(m) = costlist(id);
    
    chosen(m,:) = odemax(id,:);
    
end

%fitness function (unit if aa/h)
F = @(x) sqrt(12.5e9)*sqrt(x) - x;
%c_factor to convert aa to Ternary complex
kelongf_all = @(kmax,KM,y,Alist,b,c) T3_A(Alist,b,c).*3600.*kmax.*c_factor.*y(:,3)./(c_factor*y(:,3)+KM.*1e-6.*NA.*V(y)*1e-15);
%metabolic rate
kmetf_all = @(y,kmet) kmet./(1+(y(:,3)./(Ke*NA*V(y)*1e-15)).^2);
%input flux
JM_f = @(y,kmet) y(:,2).*kmetf_all(y,kmet);
%output flux
JS_f = @(y,kmax,KM,Alist,b,c) y(:,1).*kelongf_all(kmax,KM,y,Alist,b,c);
%degradation flux
Jd_f = @(y,Alist) degrade(Alist).*(NR*y(:,1)+NP*y(:,2));

%check flux balances
ratio_result = cell(length(kmetsample),1);
for m = 1 : length(kmetsample)
    % m = 4;
    temp = NaN*ones(length(flux_kmet{m}(:,1)),5);
    % JM - degradation flux = lambda * M
    temp(:,1) = flux_kmet{m}(:,1)./(flux_kmet{m}(:,3)+odemax_kmet{m}(:,4).*odemax_kmet{m}(:,11));
    % R*kelong + lambda*aa = JM
    temp(:,2) = flux_kmet{m}(:,1)./(flux_kmet{m}(:,2)+odemax_kmet{m}(:,4).*odemax_kmet{m}(:,3));

    %input flux check: numerical vs analytical 
    temp(:,3) = flux_kmet{m}(:,1)./JM_f(odemax_kmet{m},kmetsample(m));
    %synthesis flux check
    temp(:,4) = flux_kmet{m}(:,2)./JS_f(odemax_kmet{m},ksynmaxlist,KMefflist,Alist,blist(m),clist(m));
    %degradation flux check
    temp(:,5) = flux_kmet{m}(:,3)./Jd_f(odemax_kmet{m},Alist);
    
    ratio_result{m} = temp;
end

figure
for m = 1 : length(Alist)
    plot(odestore{m}(:,5),odestore{m}(:,4))
    hold on
end
legend('A = 1.96e4','A = 5.53e3','A = 1.90e3','A = 580','A = 92')
xlabel('Ribosome allocation')
ylabel('Growth rate (1/h)')


%nutrient quality vs growth rate vs fitness
%fitness function F = axlambda

figure
yyaxis left 
plot(flux_nutrient(:,1), chosen(:,4), '.-','MarkerFace',[0, 0.4470, 0.7410], 'MarkerSize',15)
ylabel('Growth rate (1/h)')
ylim([0 2])
yticks(0:0.4:2)
yyaxis right
% plot(nutrient, 2.5*chosen(:,4) - flux_nutrient(:,1)./2e9, 'sq-','MarkerFace',[0.8500, 0.3250, 0.0980],'MarkerSize',5)
% plot(flux_nutrient(:,1), 2.5*sqrt(flux_nutrient(:,1)/Mx) - flux_nutrient(:,1)./2e9, 'sq-','MarkerFace',[0.8500, 0.3250, 0.0980],'MarkerSize',5)
plot(flux_nutrient(:,1), F(flux_nutrient(:,1))/max(F(flux_nutrient(:,1))), 'sq-','MarkerFace',[0.8500, 0.3250, 0.0980],'MarkerSize',5)
% scatter(nutrient, chosen(:,4)./flux_nutrient(:,1), 'sq',[],'Color',[0.8500, 0.3250, 0.0980])
ylabel('Fitness')
% ylim([1 4]*1e9)
% yticks((1:1:4)*1e9)   
ylim([0 1])
yticks(0:0.2:1)  
xlabel('Intake flux (aa/h)')
xlim([0 8e9])


%plot only fitness
color_fluxlist = [77 191 237; 113 152 241; 148 114 244; 188 71 248; 255 0 255]/255; %list of color

figure
plot(flux_nutrient(:,1), F(flux_nutrient(:,1))/max(F(flux_nutrient(:,1))),'k.-')
hold on
scatter(flux_nutrient(:,1), F(flux_nutrient(:,1))/max(F(flux_nutrient(:,1))),50, color_fluxlist, 'filled')
% scatter(nutrient, chosen(:,4)./flux_nutrient(:,1), 'sq',[],'Color',[0.8500, 0.3250, 0.0980])
ylabel('Fitness')
% ylim([1 4]*1e9)
% yticks((1:1:4)*1e9)   
ylim([0 1])
yticks(0:0.2:1)  
xlabel('J_M (aa/h)')
xlim([0 8e9])
xticks((0:1:8)*1e9)
set(gca,'FontName','Helvetica','FontSize',16)

%hypothetical case where fitness max lies at same place at max growth rate
% 
% figure
% yyaxis left 
% plot(flux_nutrient(:,1), chosen(:,4), '.-','MarkerFace',[0, 0.4470, 0.7410], 'MarkerSize',15)
% ylabel('Growth rate (1/h)')
% ylim([0 2])
% yticks(0:0.4:2)
% yyaxis right
% plot(flux_nutrient(:,1), 3.3*sqrt(flux_nutrient(:,1)/Mx) - flux_nutrient(:,1)./2e9, 'sq-','MarkerFace',[0.8500, 0.3250, 0.0980],'MarkerSize',5)
% % plot((1e9:1e9:7e9), 3.3*sqrt((1e9:1e9:7e9)/Mx) - (1e9:1e9:7e9)./2e9, 'sq-','MarkerFace',[0.8500, 0.3250, 0.0980],'MarkerSize',5)
% hold on
% plot((8e9:1e9:14e9), 3.3*sqrt((8e9:1e9:14e9)/Mx) - (8e9:1e9:14e9)./2e9, 'sq-','MarkerFace',[128,128,128]/255,'MarkerSize',5,'MarkerEdgeColor', [128,128,128]/255)
% [val, id] = max(chosen(:,4));
% x_corr = [flux_nutrient(id,1) 15e9 15e9 flux_nutrient(id,1)];
% y_corr = [1 1 3.5 3.5];
% patch(x_corr,y_corr,[128,128,128]/255,'FaceAlpha',0.4,'EdgeColor','none');
% dim = [0.5 0.2 0.3 0.3];
% str = {'Physically impossible'};
% annotation('textbox',dim,'String',str,'FitBoxToText','on');
% ylabel('Fitness')
% xlabel('Intake flux (aa/h)')
% xlim([0 15e9])

figure
yyaxis left 
plot(flux_nutrient(:,1), chosen(:,4), '.-','MarkerFace',[0, 0.4470, 0.7410], 'MarkerSize',15)
ylabel('Growth rate (1/h)')
ylim([0 2])
yticks(0:0.4:2)
yyaxis right
plot(flux_nutrient(:,1), A_nutrient, 'k^-','MarkerFace','k','MarkerSize',5)
ylim([1e3 1e5])
set(gca,'YScale','log')
set(gca,'ycolor','k')
ylabel('Accuracy')
xlabel('Intake flux (aa/h)')
xlim([0 8e9])

%nutrient quality vs growth rate vs efficiency
figure
yyaxis left 
plot(nutrient, chosen(:,4), '.-','MarkerFace',[0, 0.4470, 0.7410], 'MarkerSize',15)
ylabel('Growth rate (1/h)')
yyaxis right
plot(nutrient, chosen(:,4)./flux_nutrient(:,1), 'sq-','MarkerFace',[0.8500, 0.3250, 0.0980],'MarkerSize',5)
% scatter(nutrient, chosen(:,4)./flux_nutrient(:,1), 'sq',[],'Color',[0.8500, 0.3250, 0.0980])
ylabel('Efficiecyn \lambda/J_M (1/aa)')
xlabel('Nutrient quality')

%nutrient quality vs growth rate vs flux
figure
yyaxis left 
plot(nutrient, chosen(:,4), '.-','MarkerFace',[0, 0.4470, 0.7410], 'MarkerSize',15)
ylabel('Growth rate (1/h)')
ylim([0 2])
yyaxis right
bar(nutrient, [chosen(:,4).*chosen(:,11) flux_nutrient(:,3) ], 'stacked')
hold on
scatter(nutrient, flux_nutrient(:,1),[],[0.8500, 0.3250, 0.0980],'sq','filled')
ylabel('Flux (aa/h)')
xlabel('Nutrient quality')
ylim([0 1e10])
legend('Growth rate','Growth flux','Degradation flux','Input flux')

%growth rate vs accuracy
figure
for i = 1 : length(kmetsample)
   plot(Alist, odemax_kmet{i}(:,4), '.-', 'Color', color_list(i,:)) 
   hold on
   scatter(A_nutrient(i),chosen(i,4),[],color_list(i,:),'filled')
end
set(gca,'XScale','log')
xlabel('Accuracy')
ylabel('Growth rate (1/h)')

%nutrient quality vs growth rate and accuracy 
figure
yyaxis left 
plot(nutrient, chosen(:,4), '.-','MarkerFace',[0, 0.4470, 0.7410], 'MarkerSize',15)
ylabel('Growth rate (1/h)')
yyaxis right
plot(nutrient, A_nutrient, 'sq-','MarkerFace',[0.8500, 0.3250, 0.0980],'MarkerSize',5)
xlabel('Nutrient quality')
ylabel('Accuracy')
set(gca,'YScale','log')

%plot optimal flux balance
figure
yyaxis left
bar(chosen(:,4), [chosen(:,4).*chosen(:,11) flux_nutrient(:,3) ], 'stacked')
hold on
scatter(chosen(:,4), flux_nutrient(:,1),[],'m','filled')
ylabel('Flux (aa/h)')
yyaxis right
plot(chosen(:,4), A_nutrient, 'sq-', 'MarkerFace',[0.8500, 0.3250, 0.0980] ,'MarkerSize',5)
xlabel('Growth rate (1/h)')
ylabel('Accuracy')
legend('Growth flux','Degradation flux','Input flux','Accuracy')
set(gca,'YScale','log')


figure
b = bar(chosen(:,4), [chosen(:,4).*chosen(:,11) flux_nutrient(:,3) ], 'stacked')
b(2).FaceColor = [0.9290    0.6940    0.1250];
hold on
scatter(chosen(:,4), flux_nutrient(:,1),[],'m','filled')
ylabel('Flux (aa/h)')
xlabel('Growth rate (1/h)')
legend('Growth flux','Degradation flux','Input flux')
axes('Position',[.2 .5 .3 .3])   %[xstart ystart xend-xstart yend-ystart ]
box on
bar(chosen(:,4),flux_nutrient(:,3),'FaceColor',[0.9290    0.6940    0.1250])
xlim([0 2.5])
xticks(0:0.5:2.5)
ylabel('Flux (aa/h)')
xlabel('Growth rate (1/h)')

namelist = {'R','P','aa'};
masslist = [NR NP 1];

for i = 1 : length(namelist)
    
    yL = sprintf('%s^* (molecules/cell)',namelist{i}); 
    yR = ['\phi_' sprintf('{%s}^*',namelist{i})];       
    figure
    yyaxis left
    plot( chosen(:,4),   chosen(:,i),'.-','MarkerSize',20)
    ylabel(yL)
    xlabel('Growth rate (1/h)')
    yyaxis right
    plot( chosen(:,4),  chosen(:,i+6),'sq-')
    ylabel(yR)        
end

figure
plot(chosen(:,4), kelong_nutrient/3600, '.-')
hold on
plot(chosen(:,4), koriginal_nutrient/3600, 'sq-')
xlabel('Growth rate (1/h)')
ylabel('Synthesis speed (1/s)')
legend('After recycling','No recycling')



%% compare the flux with metabolic fluxes from PNAS paper
mW = readtable('data_scale.xlsx');
mW = table2array(mW);
%1st col: cell mass (g) 
%2nd col: metabolic flux (W)
xval = (10.^(-14:0.1:-10))'; % g

color_fluxlist = [77 191 237; 113 152 241; 148 114 244; 188 71 248; 255 0 255]/255; %list of color

%convert from gram to aa 
figure
scatter(mW(:,1)*dwconv*RPpercent/maa,mW(:,2)*nATP*3600/nconv,[],[169 169 169]/255,'sq','filled')
hold on
scatter(mW(32,1)*dwconv*RPpercent/maa,mW(32,2)*nATP*3600/nconv,50,[255 153 85]/255,'filled') %Ecoli
scatter(Mpicklist, flux_nutrient(:,1),50,color_fluxlist,'filled')
set(gca,'XScale','log','YScale','log')
xlabel('M (aa)')
ylabel('J_M (aa/h)')
% legend('','E Coli')
set(gca,'FontName','Helvetica','FontSize',16)
box on

% annotation('textarrow',[0.04 0.024],[0.85 0.85] ,'String','Optimal intake flux')
box on
ticks = get(gca,'XTickLabel')
expression = '\d*\^\{(\-?\d*)\}'; % Dynamic Regular Expressions
replace = '$1'; 
exponents_X = regexprep(ticks,expression,replace);
ticks = get(gca,'YTickLabel')
exponents_Y = regexprep(ticks,expression,replace);
set(gca,'XTickLabel',exponents_X,'YTickLabel',exponents_Y)

%weighted average mass of 1 amino acid 
maa = 110*1.66e-24; % g = 110 Da

%% check time evolution of ode for max growth rate of each nutrient 
store_recheck = NaN*ones(length(kmetsample),11);
m = 5;
y0 = chosen(m,1:3);
ksynmaxT3 = ksynmax_nutrient(m)*3600;
kmet = kmetsample(m);
KMeff = KM_nutrient(m)*1e-6; %M
par = [ksynmaxT3 kmetsample(m) chosen(m,5) NR NP A_nutrient(m) blist(m) clist(m)];
Mth = chosen(m,11);
A = A_nutrient(m);
%cell mass (aa)
Mf = @(y) NR*y(1)+NP*y(2)+y(3);
%volume (um^3)
Vf = @(y) Mf(y)/rho; %um^3 
%R mass fraction
Rmfc = @(y) NR*y(1)./Mf(y);

%P mass fraction
Pmfc = @(y) NP*y(2)./Mf(y);

%aa mass fraction
aamfc = @(y) 1 - Rmfc(y) - Pmfc(y);
T3convertf_T3 = @(y) c_factor*y(3);

ktransf_T3 = @(y) ksynmaxT3*T3convertf_T3(y)./(T3convertf_T3(y)+KMeff*NA*Vf(y)*1e-15);

flist = {Mf,Vf,ktransf_T3,T3convertf_T3};

%metabolic rate
kmetf = @(y) kmet./(1+(y(3)/(Ke*NA*Vf(y)*1e-15))^2);
flist{5} = kmetf; 
flist{6} = qA; 
flist{7} = T3_A;
flist{8} = degrade;
        
k = 1;
t_tot = 0;
y_tot = y0;
tspan = [0 200];  
gr = 1;
Opt = odeset('Events', @(t,y) myEvent_growth_opt2(t,y,Mth,par));
%run to reach exponential steady state 
while k <= nSim
    [t,y, te, ye, ie] = ode15s(@(t,y)  ode_growth_linkingT3A(t,y,par,flist),tspan,y0, Opt);
    if isempty(ye) ~= 1
        t_tot = [t_tot ; [t;te] + t_tot(end)];
        y_tot = [ y_tot ; [y;ye] ];    
        y0 = ye/2;
        k = k + 1;
    else
        gr = NaN;
        break;
    end
end
%     semilogy(t_tot,y_tot)
if isnan(gr) ~= 1 
% [t,y,te,ye,ie] = ode15s(@(t,y)  ode_growth_linkingT3A(t,y,par,flist),[0 2*te],y0,Opt);
[t,y,te,ye,ie] = ode15s(@(t,y)  ode_growth_linkingT3A(t,y,par,flist),[0:1e-4:2*te],y0,Opt);
chosen(m,4)/(log(2)/t(end))
%     semilogy(t,y)
store_recheck(m,:) = [y0 log(2)/t(end) chosen(m,5) 1 - chosen(m,5) Rmfc(y0) Pmfc(y0) aamfc(y0) T3convertf_T3(y0) Mth];
end

%evaluate dy/dt 
dy_dt = NaN*ones(length(t),3);
for i = 1:length(t)
    dy_dt(i,:) = ode_growth_linkingT3A(t(i),y(i,:),par,flist);
end

%time dependent input flux 
JM = y(:,2)*kmet./(1+(y(:,3)./(Ke*V(y)*NA*1e-15)).^2);
%time dependent synthesis flux 
JS = y(:,1).*ksynmaxT3.*c_factor.*y(:,3)./(c_factor.*y(:,3)+KMeff*NA*V(y)*1e-15).*T3_A(A,b,c);

flux_nutrient(m,1)/JM(1)
flux_nutrient(m,2)/JS(1)
flux_nutrient(m,3)/(NR*y(1,1)*degrade(A)+NP*y(1,2)*degrade(A))

JM./(JS+dy_dt(:,3))
JS./(NR*dy_dt(:,1)+NP*dy_dt(:,2)+NR*y(:,1)*degrade(A)+NP*y(:,2)*degrade(A))
JM./(dy_dt(:,3)+NR*dy_dt(:,1)+NP*dy_dt(:,2)+NR*y(:,1)*degrade(A)+NP*y(:,2)*degrade(A))

hfit = fit(M(y),M(dy_dt),'poly1')
hfit.p1/((log(2)/t(end)))
flux_nutrient(m,1)./(hfit.p1*chosen(m,11)+flux_nutrient(m,3))

plot(hfit,M(y),M(dy_dt))

hfit = fit(M(y),JM-(NR*y(:,1)*degrade(A)+NP*y(:,2)*degrade(A)),'poly1');
hfit.p1/((log(2)/t(end)))
plot(hfit,M(y),JM-(NR*y(:,1)*degrade(A)+NP*y(:,2)*degrade(A)))

figure
plot(t,JM)
hold on
plot(t,M(dy_dt)+NR*y(:,1)*degrade(A)+NP*y(:,2)*degrade(A))
% plot(t,NR*dy_dt(:,1)+NP*dy_dt(:,2)+dy_dt(:,3))
legend('J_M','J_d+dM/dt')


figure
plot(t,JM)
hold on
plot(t,JS)
plot(t,NR*dy_dt(:,1)+NP*dy_dt(:,2)+dy_dt(:,3))
legend('J_M','J_S','dM/dt')
saveas(gca,fullfile(yourFolder,'time_JM'),'png')


for j = 1 : length(kmetsample)
    
     
    KMeff = KM_nutrient(j)*1e-6;
    ksynmaxT3 = ksynmax_nutrient(j)*3600;
    kmet = kmetsample(j);
    b = blist(j);
    c = clist(j);
    
    odemax = odemax_kmet{j};
    [val, id] = max(odemax(:,4));  %max growth rate
    
    par2 = [ksynmaxT3 kmet odemax(id,5) NR NP A b c Ke KMeff];

        
    T3convertf_T3 = @(y) c_factor*y(:,3);
    ktransf_T3 = @(y) ksynmaxT3*T3convertf_T3(y)./(T3convertf_T3(y)+KMeff*NA*Vf(y)*1e-15);
    flist = {Mf,Vf,ktransf_T3,T3convertf_T3};
   
    %metabolic rate
    
    kmetf = @(y) kmet./(1+(y(3)/(Ke*NA*Vf(y)*1e-15))^2);
    flist{5} = kmetf; 
    flist{6} = qA; 
    flist{7} = T3_A;
    flist{8} = degrade;


        
%     Foldername = sprintf('kmet = %.2f_timeode_T3',kmetsample(j)/3600);
    Foldername = sprintf('A_kmet = %.2f',kmetsample(j)/3600); %for T3 proportional to aa
    recheckode_T3_A(id,odemax,kmetsample,flist,par2,M,Rmfc,Pmfc,aamfc,V,kelongf_all,Foldername,Vf,T3convertf_T3)
end

