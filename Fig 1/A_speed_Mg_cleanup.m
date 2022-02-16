%use data for cognate AAA and near cognate GAA
%Johansson et al, PNAS, 2012 www.pnas.org/cgi/doi/10.1073/pnas.1116480109
%based on SI of Zhang et al, RNA, 22:896-904, 2016

%assumption
kpepnc = 0.3; %1/s
khyd = 500; %1/s GTP hydrolysis is fast 
kpepc = 7; %1/s
qc = 1;  %cognate PR rejection rate - assume constant with Mg2+ ,as long as << kpepc 
%assumption taken from 10.1016/j.molcel.2005.12.018
Rtotal = 10; %uM
%total ternary complex = 100 uM 
T3c = 2; %uM 
T3nc = 15; %uM 

%for initial selection
%Johansson et al, PNAS, 2012 www.pnas.org/cgi/doi/10.1073/pnas.1116480109
kcatKMc = [60; 117; 147; 167; 180];  % cognate AAA uM^-1s^-1
kcatKMnc = [19; 66; 139; 327; 1750]; % near cognate GAA mM^-1s^-1 
kcatKMnc = kcatKMnc/1000; %convert from mM^-1 to uM^-1

%for peptide bond
kcatKMpepnc = [3.9e-4; 2.7e-3; 9.86e-3;3.67e-2; 2.5e-1];  %unit is uM^-1s^-1
kcatKMpepc = [60; 117; 147; 167; 180]; %unit is uM^-1s^-1

%based on SI of Zhang et al, RNA, 22:896-904, 2016
Mgpep = [1.3; 2.3; 3.4; 4.6; 7.5]; %mM free Mg2+ 

%from coarse-grain scheme 
qnc = (kcatKMnc./kcatKMpepnc-1)*kpepnc; 

%concentration distribution
syms R Actc Actnc PRc PRnc

Actcsol_store = NaN*ones(length(qnc),1);
Actncsol_store = NaN*ones(length(qnc),1);
Rsol_store = NaN*ones(length(qnc),1);
PRcsol_store = NaN*ones(length(qnc),1);
PRncsol_store = NaN*ones(length(qnc),1);

R_ana = NaN*ones(length(qnc),1);
Actc_ana = NaN*ones(length(qnc),1);
Actnc_ana = NaN*ones(length(qnc),1);
PRc_ana = NaN*ones(length(qnc),1);
PRnc_ana = NaN*ones(length(qnc),1);

R_approx = NaN*ones(length(qnc),1);
Actc_approx = NaN*ones(length(qnc),1);
Actnc_approx = NaN*ones(length(qnc),1);
PRc_approx = NaN*ones(length(qnc),1);
PRnc_approx = NaN*ones(length(qnc),1);

for i = 1 : length(kcatKMc)
    %equation at steady state 
    eq1 = R * kcatKMc(i) * T3c - (khyd)*Actc == 0;
    eq2 = Actc * khyd - (qc + kpepc)*PRc == 0;

    eq3 =  R * kcatKMnc(i) * T3nc - (khyd)*Actnc == 0;
    eq4 = Actnc * khyd - (qnc(i) + kpepnc)*PRnc== 0;

    %eq5 =  (qc + kpepc)*PRc + (qnc + kpepnc)*PRnc- R * (k11 + k12) == 0

    eq7 = R + Actc + Actnc + PRc + PRnc == Rtotal ;

    eqn = [eq1, eq2, eq3, eq4, eq7];

    [Rsol, Actcsol, Actncsol, PRcsol, PRncsol ] = solve(eqn, [R Actc Actnc PRc PRnc]);
    Rsol_store(i,1) = Rsol;
    Actcsol_store(i,1) = Actcsol;
    Actncsol_store(i,1) = Actncsol;
    PRcsol_store(i,1) = PRcsol; 
    PRncsol_store(i,1) = PRncsol;
    
    %analytical solution
    R_ana(i,1) = Rtotal/(1+T3c*kcatKMc(i)/khyd*(1+khyd/(qc+kpepc))+T3nc*kcatKMnc(i)/khyd*(1+khyd/(qnc(i)+kpepnc)));
    Actc_ana(i,1) = Rtotal*(T3c*kcatKMc(i)/khyd)/(1+T3c*kcatKMc(i)/khyd*(1+khyd/(qc+kpepc))+T3nc*kcatKMnc(i)/khyd*(1+khyd/(qnc(i)+kpepnc)));
    Actnc_ana(i,1) = Rtotal*(T3nc*kcatKMnc(i)/khyd)/(1+T3c*kcatKMc(i)/khyd*(1+khyd/(qc+kpepc))+T3nc*kcatKMnc(i)/khyd*(1+khyd/(qnc(i)+kpepnc)));
    PRc_ana(i,1) = Rtotal*(T3c*kcatKMc(i)/khyd*khyd/(qc+kpepc))/(1+T3c*kcatKMc(i)/khyd*(1+khyd/(qc+kpepc))+T3nc*kcatKMnc(i)/khyd*(1+khyd/(qnc(i)+kpepnc)));
    PRnc_ana(i,1) = Rtotal*(T3nc*kcatKMnc(i)/khyd*khyd/(qnc(i)+kpepnc))/(1+T3c*kcatKMc(i)/khyd*(1+khyd/(qc+kpepc))+T3nc*kcatKMnc(i)/khyd*(1+khyd/(qnc(i)+kpepnc)));
    
    %approximation
    Zapprox = 1+T3c*kcatKMc(i)/kpepc+T3nc*kcatKMnc(i)/qnc(i);
    R_approx(i,1) = Rtotal/Zapprox;
    Actc_approx(i,1) = Rtotal*(T3c*kcatKMc(i)/khyd)/Zapprox;
    Actnc_approx(i,1) = Rtotal*(T3nc*kcatKMnc(i)/khyd)/Zapprox;
    PRc_approx(i,1) = Rtotal*(T3c*kcatKMc(i)/kpepc)/Zapprox;
    PRnc_approx(i,1) = Rtotal*(T3nc*kcatKMnc(i)/qnc(i))/Zapprox;


end

%compare analytical and numerical

plot(Mgpep, Rsol_store, Mgpep, R_ana,'*')
plot(Mgpep, Actcsol_store, Mgpep, Actc_ana,'*')
plot(Mgpep, Actncsol_store, Mgpep, Actnc_ana,'*')
plot(Mgpep, PRcsol_store, Mgpep, PRc_ana,'*')
plot(Mgpep, PRncsol_store, Mgpep, PRnc_ana,'*')


scatter(Mgpep, PRc_ana, 'filled')
hold on
plot(Mgpep, PRc_approx,'*-')
xlabel('Free [Mg^{2+}] (mM)')
ylabel('[PR^{c}] (\muM)')
legend('Exact','Approximation')
xlim([1 8])

scatter(Mgpep, PRnc_ana,'filled')
hold on
plot(Mgpep, PRnc_approx,'*-')
xlabel('Free [Mg^{2+}] (mM)')
ylabel('[PR^{nc}] (\muM)')
legend('Exact','Approximation')
xlim([1 8])

%justify rate approximation to get ksynmax and KMeff
ksyneff_ana = (PRc_ana*kpepc+PRnc_ana*kpepnc)/Rtotal;
ksyneff_approx1 =  PRc_approx*kpepc/Rtotal;
ksyneff_approx2 =  (PRc_approx*kpepc+PRnc_approx*kpepnc)/Rtotal;

figure
scatter(Mgpep, ksyneff_ana,'filled')
hold on
plot(Mgpep, ksyneff_approx1,'*-')
% plot(Mgpep, ksyneff_approx2,'o-')
xlabel('Free [Mg^{2+}] (mM)')
ylabel('k^{eff}_{syn} (aa/(R.s)')
xlim([1 8])

Jtotal_normalized = (PRcsol_store*kpepc + PRncsol_store*kpepnc)/Rtotal;
Atotal = PRcsol_store*kpepc./(PRncsol_store*kpepnc);
JIS_normalized =  (Actcsol_store + Actncsol_store)*khyd/Rtotal;
Q_normalized = (PRcsol_store*qc + PRncsol_store.*qnc)/Rtotal;
AIS = Actcsol_store./Actncsol_store;
Aproof = (kpepc/(kpepc+qc))./(kpepnc./(kpepnc+qnc));
Jtotalc_normalized = PRcsol_store*kpepc/Rtotal;
Jtotalnc_normalized = PRncsol_store*kpepnc/Rtotal;
Qc_normalized = PRcsol_store*qc/Rtotal;
Qnc_normalized = PRncsol_store.*qnc/Rtotal;
JISc_normalized =  Actcsol_store*khyd/Rtotal;
JISnc_normalized =  Actncsol_store*khyd/Rtotal;
Jactc_normalized =  Rsol_store.*kcatKMc*T3c/Rtotal;
Jactnc_normalized = Rsol_store.*kcatKMnc*T3nc/Rtotal;

%ksyneff expressed in Michaelis - Menten form 
fc = 0.02;
fnc = 0.15;
T3 = 100; %uM 
ksynmax = @(kfc,kpepc,kfnc,krejnc) fc*kfc./(fc*kfc/kpepc+fnc*kfnc./krejnc);
KMeff = @(kfc,kpepc,kfnc,krejnc) 1./(fc*kfc/kpepc+fnc*kfnc./krejnc);
Accuracy = @(kfc,kpepc,kfnc,krejnc,krejc) fc*kfc./(fnc*kfnc).*(kpepc/(krejc+kpepc))./(kpepnc./(krejnc+kpepnc));
plot(Mgpep, ksynmax(kcatKMc,kpepc,kcatKMnc,qnc))
plot(Mgpep, KMeff(kcatKMc,kpepc,kcatKMnc,qnc))

plot(Mgpep, ksynmax(kcatKMc,kpepc,kcatKMnc,qnc)*T3./(T3+KMeff(kcatKMc,kpepc,kcatKMnc,qnc)))
hold on
scatter(Mgpep, Jtotal_normalized)

plot(Mgpep, Accuracy(kcatKMc,kpepc,kcatKMnc,qnc,qc))
hold on
scatter(Mgpep, Atotal)

%illustrate effect of competition - put in paper SI 
i = 3;
kfc = kcatKMc(i);
kfnc = kcatKMnc(i);
krejnc = qnc(i);
fc = 0.02; 
T3 = 100; %uM 
ksm = @(fnc) fc*kfc./(fc*kfc/kpepc+fnc*kfnc./krejnc);
Km =  @(fnc) 1./(fc*kfc/kpepc+fnc*kfnc./krejnc);

fnclist = (0:0.01:0.5)';

figure
yyaxis left
plot(fnclist/fc, ksm(fnclist))
ylabel('k^{max}_{syn} (s^{-1})')
yyaxis right
plot(fnclist/fc, Km(fnclist))
ylabel('K^{eff}_{M} (\muM)')
xlabel('Near cognates/cognates ratio')

figure
plot(fnclist/fc, ksm(fnclist)*T3./(T3+Km(fnclist)))
xlabel('Near cognates/cognates ratio')
ylabel('k^{eff}_{syn} (s^{-1})')


%approach 1 - rate approximation with Mg2+ 
[fit_kcatKMc,fit_kcatKMnc,fb,fa] = appproach1_rateapprox(Mgpep, kcatKMc, kcatKMnc, qnc, kpepnc, qc, kpepc, khyd)
Mgseries = 0:0.1:8;

%plot rejection rate fit
figure
scatter(Mgpep, qnc, 'filled')
hold on
plot(Mgseries, fit_kcatKMnc(Mgseries)./fb(Mgseries)-kpepnc)
xlabel('Free [Mg^{2+}] (mM)')
ylabel('k^{rej}_{nc} (aa/(R.s)')
xlim([1 8])

PRc_approx = @(x,T3c, T3nc) T3c*fa(x)./(1+T3c*fa(x)+T3nc*fb(x));
PRnc_approx = @(x,T3c, T3nc) T3nc*fb(x)./(1+T3c*fa(x)+T3nc*fb(x));
R_approx = @(x,T3c, T3nc) 1./(1+T3c*fa(x)+T3nc*fb(x));

%approach 2 - concentration approximation with Mg2+ 
[ncog_fit,r_fit, cog_fit]= approach2_Capprox(Mgpep,Rtotal,Actcsol_store, Actncsol_store,Rsol_store,PRcsol_store, PRncsol_store)

%some plots to check quality of fitting 

figure
scatter(Mgpep, Rsol_store*100/Rtotal,'filled')
hold on
plot(Mgseries, Rtotal*R_approx(Mgseries,T3c, T3nc)*100/Rtotal,'k')
plot(Mgseries, r_fit(Mgseries)*100/Rtotal,'k--')
xlabel('Free [Mg^{2+}] (mM)')
ylabel('[R]_{standby} fraction (%)')
legend('Model','Rate approximation','Concentration approximation','Location','NorthEast')
xlim([1 8])
% saveas(gca,fullfile(yourFolder,'kcatKM'),'png')
saveas(gca,'Rstandby_fittingcompare','png')

figure
scatter(Mgpep, PRcsol_store*100/Rtotal,'filled')
hold on
plot(Mgseries, Rtotal*PRc_approx(Mgseries,T3c, T3nc)*100/Rtotal,'k')
plot(Mgseries, cog_fit(Rtotal,Mgseries)*100/Rtotal,'k--')
xlabel('Free [Mg^{2+}] (mM)')
ylabel('[PR]_{c} fraction (%)')
legend('Model','Rate approximation','Concentration approximation','Location','SouthEast')
xlim([1 8])
% saveas(gca,fullfile(yourFolder,'kcatKM'),'png')
saveas(gca,'PRc_fittingcompare','png')

figure
scatter(Mgpep, PRncsol_store*100/Rtotal,'filled')
hold on
plot(Mgseries, Rtotal*PRnc_approx(Mgseries,T3c, T3nc)*100/Rtotal,'k')
plot(Mgseries, ncog_fit(Mgseries)*100/Rtotal,'k--')
xlabel('Free [Mg^{2+}] (mM)')
ylabel('[PR]_{nc} fraction (%)')
legend('Model','Rate approximation','Concentration approximation','Location','NorthWest')
xlim([1 8])
saveas(gca,'PRnc_fittingcompare','png')

%rate approximation
figure
plot(Mgseries, 100*PRc_approx(Mgseries,T3c, T3nc))
hold on
plot(Mgseries, 100*(PRc_approx(Mgseries,T3c, T3nc)+PRnc_approx(Mgseries,T3c, T3nc)))
ylabel('Ribosome partition (%)')
xlabel('Free [Mg^{2+}] (mM)')
ylim([70 100])
xlim([1 8])

%concentration approximation
figure
plot(Mgseries, cog_fit(Rtotal,Mgseries)/Rtotal*100)
hold on
plot(Mgseries, (cog_fit(Rtotal,Mgseries)+ncog_fit(Mgseries))/Rtotal*100)
ylabel('Ribosome partition (%)')
xlabel('Free [Mg^{2+}] (mM)')
ylim([70 100])
xlim([1 8])

%run A_Mg_testerror.m and function therein 
run('A_Mg_testerror.m')
close all
figure
% scatter(Mgpep, Jtotal_normalized,'filled')
errorbar(Mgpep, Jtotal_normalized, derry, uerry,'s','MarkerFaceColor',[0 0.45 0.74])
hold on
plot(Mgseries, PRc_approx(Mgseries,T3c, T3nc)*kpepc + PRnc_approx(Mgseries,T3c, T3nc)*kpepnc,'k')
plot(Mgseries, (cog_fit(Rtotal,Mgseries)*kpepc + ncog_fit(Mgseries)*kpepnc)/Rtotal,'k--')
xlabel('Free [Mg^{2+}] (mM)')
ylabel('k_{syn}^{eff} (s^{-1})')
legend('Model','Rate parameterization','Concentration parameterization','Location','SouthWest')
xlim([1 8])
ylim([5 7])
saveas(gca,'ksyneff_fittingcompare','png')
saveas(gca,'ksyneff_fittingcompare.eps','eps')

close all



%% plot figure for publication
%approximation flux 
J_1 = Rtotal*PRc_approx(Mgseries,T3c, T3nc)*kpepc + Rtotal*PRnc_approx(Mgseries,T3c, T3nc)*kpepnc;
J_2 = ncog_fit(Mgseries)*kpepnc + cog_fit(Rtotal,Mgseries)*kpepc;

%normalized flux (aka effective synthesis speed)
Japprox_1 = J_1/Rtotal;
Japprox_2 = J_2/Rtotal;
%accuracy
A_1 = (PRc_approx(Mgseries,T3c, T3nc)*kpepc)./(PRnc_approx(Mgseries,T3c, T3nc)*kpepnc);
A_2 = (cog_fit(Rtotal,Mgseries)*kpepc)./(ncog_fit(Mgseries)*kpepnc);

%simplified rate-approximation normalized peptide flux further 
a = 23.11;
b = 0.13;
c = 1.11;
p = 7.19e-4;
q = 1;      
T3 = 100; %T3 unit is uM 
fc = 0.02; %cognate fraction
fnc = 0.15; %near cognate fraction
% J_simpified = @(x,fc,fnc,T3) T3*fc*a./(b+exp(-c*x)+T3*(fc*a/kpepc+fnc*p*(b*exp(q*x)+exp((q-c)*x)))); 
J_simpified = @(x,fc,fnc,T3) T3*fc*a./(b+exp(-c*x)+T3*(fc*a/kpepc+fnc*p*(b*exp(q*x)))); 
Mgpeak = 1/(q+c)*log(c/(T3*fnc*p*b*q));

%figure in SI to compare ksyneff between two parameterization 
figure
scatter(Mgpep, Jtotal_normalized,'filled')
hold on
plot(Mgseries,J_1/Rtotal,'k',Mgseries, J_2/Rtotal,'k--')
% plot(Mgseries,J_simpified(Mgseries,0.02,0.15,100),'k--')
% plot([Mgpeak Mgpeak],[3 7]) 
xlabel('Free [Mg^{2+}] (mM)')
ylabel('k^{eff}_{syn} (s^{-1})')
legend('Model','Rate parameterization','Concentration parameterization')
xlim([1 8])

figure
plot(Mgseries,A_1,Mgseries, A_2)
xlabel('Free [Mg^{2+}] (mM)')
ylabel('Accuracy')
legend('Rate approximation','Concentration approximation')
xlim([1 8])


%plot Rstandby, PRc and PRnc distribution 
con_fit = r_fit(Mgseries)+cog_fit(Rtotal,Mgseries)+ncog_fit(Mgseries);

figure
area(Mgseries,con_fit*100/Rtotal,'FaceColor','w','EdgeColor','k','FaceAlpha',0.3)
hold on
area(Mgseries,cog_fit(Rtotal,Mgseries)*100/Rtotal,'FaceColor',[5 115 86]/255,'EdgeColor',[5 115 86]/255,'FaceAlpha',0.3)
p = fill([Mgseries, fliplr(Mgseries)], [cog_fit(Rtotal,Mgseries), fliplr((cog_fit(Rtotal,Mgseries)+ncog_fit(Mgseries)))]*100/Rtotal, [139,0,0]/255)
p.EdgeColor = [139,0,0]/255;
p.FaceAlpha = 0.3;
% ylim([6.5 10])
ylim([70 100])
xlim([1 8])
xlabel('Free [Mg^{2+}] (mM)')
ylabel('Concentration (\muM)')
legend('R_{standby}','PR^c','PR^{nc}','Location','SouthWest')
% saveas(gca,fullfile(yourFolder,sprintf('Mg_Rconcentrationapprox %d',i)),'png')


%3D plot 
[val, idx] = max(J_2);
pick = [10 idx 50];

g = figure

plot3(Mgseries,Japprox_2,A_2,'LineWidth',5,'Color',[50 50 50]/255)
hold on
xlim1 = 0.5; xlim2 = 8;
ylim1 = 5.2; ylim2 = 7;
zlim1 = -250; zlim2 = 6000;
xlim([xlim1 xlim2])
ylim([ylim1 ylim2])
zlim([zlim1 zlim2])

scatter3(Mgseries(idx),Japprox_2(idx),A_2(idx),'o','Filled','MarkerFaceColor',[50 50 50]/255,'MarkerEdgeColor',[50 50 50]/255,'LineWidth',5)
plot3(Mgseries,Japprox_2,zlim1*ones(length(Mgseries),1),'m','LineWidth',1);
line([Mgseries(idx) Mgseries(idx)], [val ylim1], [zlim1 zlim1],'LineStyle','--','Color','m');

scatter3(Mgseries(idx),Japprox_2(idx),zlim1,'mo','filled')
scatter3(xlim1,Japprox_2(idx),A_2(idx),'o','filled','MarkerFaceColor',[0/255 0/255 205/255],'MarkerEdgeColor',[0/255 0/255 205/255])
scatter3(Mgseries(idx),ylim1,A_2(idx),'o','filled','MarkerFaceColor',[153/255 51/255 0/255],'MarkerEdgeColor',[153/255 51/255 0/255])

line(xlim1*ones(length(Mgseries),1),Japprox_2,A_2,'LineStyle','-','Color',[0/255 0/255 205/255],'LineWidth',1);

line([xlim1 xlim1],[ylim1 Japprox_2(idx)],[A_2(idx) A_2(idx)],'LineStyle','--','Color',[128/255 128/255 128/255]);
line(Mgseries,ylim1*ones(length(Mgseries),1),A_2,'LineStyle','-','Color',[153/255 51/255 0/255],'LineWidth',1);

for i = 1 : 3    
    line([Mgseries(pick(i)) Mgseries(pick(i))],[Japprox_2(pick(i)) Japprox_2(pick(i))],[A_2(pick(i)) 10],'LineStyle','--','Color','m')
    line([Mgseries(pick(i)) Mgseries(pick(i))],[ylim1 Japprox_2(pick(i))],[A_2(pick(i)) A_2(pick(i))],'LineStyle','--','Color',[128/255 128/255 128/255])
    line([xlim1 Mgseries(pick(i))],[Japprox_2(pick(i)) Japprox_2(pick(i))],[A_2(pick(i)) A_2(pick(i))],'LineStyle','--','Color',[0 0/255 205/255])
end
view([-118 33]);
set(gca,'XDir','Reverse')
set(gca,'XTickLabel',[2 4 6 8])
box on
grid on
set(gca,'ZMinorGrid','off','ZMinorTick','off')
set(gca,'FontName','Arial')
set(gca,'DataAspectratio',[4 1 3500])
%set(gca,'XTickLabel',1:8);
%set(gca,'DataAspectratio',[4 1 6660]) for log
xlabel('[Mg^{2+}] (mM)')
ylabel('k_{syn}^{eff} (aa/sec)')
zlabel('Accuracy')
