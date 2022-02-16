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

%try different set 
kcatKMc = [40; 60; 107; 135; 148; 170];  % cognate AAA uM^-1s^-1
kcatKMnc = [5; 19; 66; 119; 276; 1450]; % near cognate GAA mM^-1s^-1 
kcatKMnc = kcatKMnc/1000; %convert from mM^-1 to uM^-1

kcatKMpepnc = [1e-4; 3.9e-4; 2.4e-3; 4.86e-3;3.67e-2; 2.5e-1];  %unit is uM^-1s^-1
kcatKMpepc = [40; 60; 117; 147; 167; 180]; %unit is uM^-1s^-1

qnc = (kcatKMnc./kcatKMpepnc-1)*kpepnc; 

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

ksynmaxlist = ksynmax(kcatKMc,kpepc,kcatKMnc,qnc);
KMefflist = KMeff(kcatKMc,kpepc,kcatKMnc,qnc);
Alist = Accuracy(kcatKMc,kpepc,kcatKMnc,qnc,qc);
grlist = lambda(ksynmaxlist*T3./(T3+KMefflist)*3600);
costlist = cost_f(kcatKMc,kpepc,kcatKMnc,qnc,qc);

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


%% relate competition effect on growth rate

% remember to run bionumers.m

run('bionumers.m')

close all

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

%% explicit simulation instead of putting in function

Mpick = 0.8476e9; %aa
nSim = 100;
uRlist = 0.01:0.01:0.9;
% Kaselect = 2.5e-4*KMefflist/KMefflist(3);
Kaselect = 2.5e-4*ones(length(ksynmaxlist),1); 

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

color_list = [77 191 237; 113 152 241; 148 114 244; 188 71 248; 255 0 255]/255; %list of color

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


       T3convertf_T3 = @(y) c_factor*y(3);
      ktransf_T3 = @(y) ksynmaxT3*T3convertf_T3(y)./(T3convertf_T3(y)+KMeff*NA*Vf(y)*1e-15);

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

            y0 = [R0 P0 aa0];
                par = [ksynmaxT3 kmet uRlist(i) NR NP A b c];

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

    %analytical growth rate
    aa_ana = @(ksynmax,KMeff,A,b,c) (NR*(KMeff*1e-6/c_factor)*kmet*Ke^2./(2*ksynmax.*qA(A).*T3_A(A,b,c)*NP)).^(1/3);
    gr_ana = @(ksynmax,KMeff,A,b,c)  aa_ana(ksynmax,KMeff,A,b,c)./(3*NR*(KMeff*1e-6/c_factor)./(2*ksynmax.*qA(A).*T3_A(A,b,c))+aa_ana(ksynmax,KMeff,A,b,c).*(NP/kmet+NR./(ksynmax.*qA(A).*T3_A(A,b,c))));

    lambda_ana = gr_ana(ksynmaxlist*3600, KMefflist, Alist, b, c);


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

%figure 2d
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

%plot optimal flux balance
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


