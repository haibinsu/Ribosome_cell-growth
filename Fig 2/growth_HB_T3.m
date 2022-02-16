%assume a relation between fastest growth rate and initial mass
%use KMeff from part 1 for ternary complex instead of rescaling 

%starting point is growth_HB.m 

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
S = R/5;
ksynmaxT3 = 22*3600; %1/h
% KM = 20e-6; %M
KM = 10e-6; %M 
KMeff = KM;

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

%mapping aa to T3
%exact mapping 
% T3convert_T3 = @(y)
% k2*S/R*KM*NA*V(y).*1e-15.*y(:,3)./(ksynmaxT3*Ka*NA*V(y)*1e-15+(ksynmaxT3-k2*S/R)*y(:,3));
%linear mapping
T3convert_T3 = @(y) k2*S/R*KM*y(:,3)./(ksynmaxT3*Ka);

%translation rate using T3
ktrans_T3 = @(y) ksynmaxT3*T3convert_T3(y)./(T3convert_T3(y)+KMeff*NA*V(y)*1e-15);

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

%mapping aa to T3
% T3convertf_T3 = @(y) k2*S/R*KM*NA*Vf(y)*1e-15*y(3)./(ksynmax*Ka*NA*Vf(y)*1e-15+(ksynmax-k2*S/R)*y(3));
T3convertf_T3 = @(y) k2*S/R*KM*y(3)./(ksynmax*Ka);
c_factor = (k2*S/R*KM)/(ksynmax*Ka)
%translation rate based on T3
ktransf_T3 = @(y) ksynmax*T3convertf_T3(y)./(T3convertf_T3(y)+KMeff*NA*Vf(y)*1e-15);

flist = {Mf,Vf,ktransf_T3,T3convertf_T3};

%% for each kmet, pick the fastest growth rate
% idchoose = [2;4;7;11;17];
% load('kmetsample')
% load('grsample')

% kmetsample = kmetsample + 0.025*3600;

%data from SI of  10.1038/nmicrobiol.2016.231 - Hwa data
load('dataextract')
%1st - 3rd col: #R/cell - #tRNA/cell - growth rate (1/h)
%4th col: peptide synthesis speed (aa/(R*s)) estimated from data_Dai.mat 
%5th col: R + P dry mass (pg) estimated from Liudata.mat 

load('d2')
%1st - 2nd col: growth rate - R mass fraction 

mW = readtable('data_scale.xlsx');
mW = table2array(mW);
%1st col: cell mass (g) 
%2nd col: metabolic flux (W)
xval = (10.^(-14:0.1:-10))'; % g

%remember to run bionumers.m

nATP = 1e19; % 1 J = 1e19 ATP molecules
%http://book.bionumbers.org/how-much-energy-is-released-in-atp-hydrolysis/
%https://www.ncbi.nlm.nih.gov/books/NBK7919/
nconv = 6; %6 ATP needed for aa incorporation into peptide chain
%www.pnas.org/cgi/doi/10.1073/pnas.1421138111 

mW_aa = NaN*ones(length(mW),2);
mW_aa(:,1) = mW(:,1)*dwconv*RPpercent/maa; %covnert from gram to aa 
mW_aa(:,2) = mW(:,2)*nATP*3600/nconv; % convert from Watt to aa/h

%old values when mistaking NP = 400 for NP = 300 in threshold function 
% idchoose = [2;7;11;17;24]; %values choose 
% load('kmetsample_old')
% load('grsample_old')
% My = 8e7;
% Mx = (2e9-8e7)/2;

Mx = 2.95e9; %from fitting data of Liu R+P vs growth rate 
My = 1.6e8;
Mx = Mx*RPpercent; % use a bit smaller so the calculated cell mass = R+P of experiment
My = My*RPpercent; % use a bit smaller so the calculated cell mass = R+P of experiment
M0f = @(x) My + Mx*x;

%for Ka = 1 mM and use [T3]/(KMeff + [T3]) = S/R*k2*[aa]/(K_a+[aa])
Mx = 0.8e9;
My = 2e8;
M0f = @(x) My + Mx*x;

%for Ka = 0.25 mM and [T3] = 2*k2*S/R*KM*[aa]./(ksynmax*Ka);
Mx = 1.65e9;
My = 2e8;
M0f = @(x) My + Mx*x;

% nstore = NaN*ones(length(kmetsample),6);
% %1 - 6th col: #R - #P - #aa - gr1 - gr2 - M0
% 
% for j = 1 : length(kmetsample)
% %     nstore(j,6) = M0f(grmaxstore(idchoose(j),1));
%     nstore(j,6) = M0f(grsample(j,1));
% 
% %     nstore(j,1:3) = maxstore3(j,1:3)*1e-6*NA*nstore(j,6)/rho*1e-15;
% %     nstore(j,4:5) = maxstore3(j,4:5);   
% end


% Mpick = [1.0907e9 1.4075e9 1.7476e9 1.91334e9 2.378e9]

%finally use this one 
% Mpick = [1.0907e9 1.4075e9 1.7476e9 1.91334e9 2.178e9];
% 
% Mpick = [0.97e9 1.24075e9 1.7476e9 1.91334e9 2.478e9];
% 
% Mpick = [1.14075e9 1.52e9 1.7476e9 1.91334e9 2.478e9];
% kmetsample = [0.07 0.1 0.12 0.14 0.2]'*3600;
% 
% Mpick = [1.04075e9 1.22e9 1.4476e9 1.691334e9 1.8478e9];
% kmetsample = [0.09 0.1 0.12 0.14 0.18]'*3600-0.02*3600;
% 
% Mpick = [1.04075e9 1.22e9 1.4476e9 1.691334e9 2.6478e9];
% kmetsample = [0.09 0.1 0.12 0.14 0.26]'*3600-0.02*3600;
% 
% Mpick = [1.04075e9 1.4476e9 1.84e9  2.2691334e9 2.6478e9];
% kmetsample = [0.09 0.12 0.18 0.22 0.26]'*3600-0.02*3600;
% 
% Mpick = [1.04075e9 1.4476e9 1.84e9 2.08e9 2.2691334e9];
% kmetsample = [0.09 0.12 0.16 0.18 0.2]'*3600-0.02*3600;
% 
% Mpick = [1.14075e9 1.4476e9 1.84e9 2.2691334e9 2.6478e9];
% kmetsample = [0.07 0.1 0.14 0.17 0.24]'*3600;

%use this one 
Mpick = [1.14075e9 1.6476e9 1.84e9 2.2691334e9 2.6478e9];
kmetsample = [0.07 0.1 0.14 0.17 0.24]'*3600;
nstore = NaN*ones(length(kmetsample),6);
%1 - 6th col: #R - #P - #aa - gr1 - gr2 - M0


%to see the effect of kmet and mass
%fix mass, vary kmet
% Mpick = [1.14075e9 1.6476e9 1.84e9 2.2691334e9 2.6478e9];
% kmetsample = [0.07 0.1 0.14 0.17 0.24]'*3600 + 0.05*3600;

%vary mass, fix kmet
% Mpick = [1.14075e9 1.6476e9 1.84e9 2.2691334e9 2.6478e9] - 0.2e9;
% kmetsample = [0.07 0.1 0.14 0.17 0.24]'*3600;


nstore(:,6) = Mpick';

nSim = 50;
uRlist = 0.01:0.01:0.9;
odestore = cell(length(kmetsample),1);
odemax = NaN*ones(length(kmetsample),11);
%1st - 6th col: R - P - aa - gr - uR - uP 
%7th - 11th col: phiR - phiP - phiaa - ternary complex T3 - total cell mass M0 

for j = 1 : length(kmetsample)
    kmet = kmetsample(j);
    %metabolic rate
    kmetf = @(y) kmet./(1+(y(3)/(Ke*NA*Vf(y)*1e-15))^2);
    flist{5} = kmetf; 
    store = NaN*ones(length(uRlist),10);
    %1st - 10th col: R - P - aa - gr - uR - uP - phiR - phiP - phiaa - T3
    
    Mth = nstore(j,6); % note that 2*Mth = mass threshold for cell division
%     Mth = 2.2048e+09; %fixed mass threshold
    parfor i =  1 : length(uRlist)

        R0 = 2e4;
        P0 = 1e5;
        M0 = 3e8; %aa
        aa0 = M0 - (NR*R0 + P0*NP);
%         R0 = nstore(j,1);
%         P0 = nstore(j,2);
%         aa0 = nstore(j,3);
        y0 = [R0 P0 aa0];
        par = [ksynmax kmet uRlist(i) NR NP];

        k = 1;
        t_tot = 0;
        y_tot = y0;
        tspan = [0 200];  
        gr = 1;
        Opt = odeset('Events', @(t,y) myEvent_growth_opt2(t,y,Mth,par));
        %run to reach exponential steady state 
        while k <= nSim
            [t,y, te, ye, ie] = ode15s(@(t,y)  ode_growth_HB_T3(t,y,par,flist),tspan,y0, Opt);
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
        [t,y,te,ye,ie] = ode15s(@(t,y)  ode_growth_HB_T3(t,y,par,flist),tspan,y0,Opt);
    %     semilogy(t,y)
        store(i,:) = [y0 log(2)/t(end) uRlist(i) 1 - uRlist(i) Rmfc(y0) Pmfc(y0) aamfc(y0) T3convertf_T3(y0)];
        end
    end
    odestore{j} = store;
    [val, idx] = max(store(:,4));
    odemax(j,:) = [store(idx,:) M(store(idx,:))];
end
   
%consistent check
%mass 
odemax(:,11)./nstore(:,6)
odemax(:,11)./Mpick'

%growth rate
% odemax(:,4)./grsample

% %analytical expression for growth rate
% ksynmax_change = ksynmax; 
% Ka_change = KMeff/(2*k2*S/R*KMeff/(ksynmax*Ka));
% %x is kmeteff
% aafr_theoryapprox = @(x) NA/(rho*1e15)*(NR*Ka_change*x*Ke^2/(2*ksynmax_change*NP)).^(1/3); %x is kmeteff
% gr_theory = @(x) (aafr_theoryapprox(x)*rho*1e15/NA)./(3*NR*Ka_change/(2*ksynmax_change)+(aafr_theoryapprox(x)*rho*1e15/NA).*(NP./x+NR/ksynmax_change));
% phiR_theory = @(x) NR*NA/(rho*1e15)*gr_theory(x).*(rho*1e15/NA-(aafr_theoryapprox(x)*rho*1e15/NA))./(ksynmax_change*(aafr_theoryapprox(x)*rho*1e15/NA)./(Ka_change+(aafr_theoryapprox(x)*rho*1e15/NA)));
% 
% scatter(kmetsample,odemax(:,4))
% hold on
% plot(0:0.1:0.3*3600,gr_theory(0:0.1:0.3*3600))
% 
% scatter(odemax(:,4),odemax(:,7))
% hold on
% plot(gr_theory(0:0.1:0.3*3600),phiR_theory(0:0.1:0.3*3600))
% 
% M0f = @(x) 1.56e9*x + 1.28e8
% figure
% scatter(gr',Mpercell*0.6)  
% hold on
% scatter(odemax(:,4),odemax(:,11),'sq')  
% plot(0:0.1:2,M0f(0:0.1:2))
% 
% figure
% scatter(odemax(:,4),odemax(:,1))
% hold on
% plot(gr_theory(0:0.1:0.3*3600),1/NR*phiR_theory(0:0.1:0.3*3600).*M0f(gr_theory(0:0.1:0.3*3600)))
% scatter(gr,Rpercell)
% scatter(dataextract(:,3),dataextract(:,1),50,'sq','MarkerEdgeColor','k','LineWidth',1.25,'MarkerEdgeAlpha',1) %[0 0.5 0] [0.49 0.18 0.56]
% 
% figure
% plot((0:0.1:0.3*3600)/3600,gr_theory(0:0.1:0.3*3600))
% hold on
% scatter(kmetsample/3600, odemax(:,4))

%T3 ratio
T3convert_T3(odemax)./odemax(:,10)

%mass fraction
mfraction = NaN*ones(length(kmetsample),3);
mfraction(:,1) = Rmf(odemax);
mfraction(:,2) = Pmf(odemax);
mfraction(:,3) = aamf(odemax);

mfraction(:,1)./odemax(:,7)

%in Bremmer and Dennis  10.1128/ecosal.5.2.3
gr = [0.6 1 1.5 2 2.5 3]*log(2); %growth rate 1/h
RNApercell = [23 44 76 128 180 214]; %ug/10^9 cells
Rpercell = [8 15 26 44 61 73]*1e3; 
kselect = [13 18 21 22 22 22]*3600; %aa/(R.h)

massfit = fit(vertcat(odemax(:,4), gr'),vertcat(odemax(:,11), Mpercell*0.6),'poly1')
Mx = massfit.p1;
My = massfit.p2;
M0f = @(x) My + Mx*x;

fit(gr',Mpercell*0.6,'poly1')  
1.62e9*odemax(:,4)+1.28e8

1.62e9*1+1.28e8

%to illustrate input fluxes drive growth
para = [Ka Ke NA ksynmax];
% Foldername = 'allkmetsameplot_T3';
Foldername = 'allkmetsameplot_T3_2'; %for T3 proportional to aa
fluxdrive_allkmetT3(odestore,para,V,M,kmetsample,odemax,nstore,Foldername,ktrans_T3,glawdata)

%in Bremmer and Dennis  10.1128/ecosal.5.2.3
gr = [0.6 1 1.5 2 2.5 3]*log(2); %growth rate 1/h
RNApercell = [23 44 76 128 180 214]; %ug/10^9 cells
Rpercell = [8 15 26 44 61 73]*1e3; 
kselect = [13 18 21 22 22 22]*3600; %aa/(R.h)

% ft = fittype( 'a*x^2+b*x+c', 'independent', 'x', 'dependent', 'y' );
% opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
% opts.Display = 'Off';
% opts.Lower = [1e3 100 1000];
% opts.Robust = 'Bisquare';
% opts.StartPoint = [1e4 5000 300];
% [gfit,gof] = fit(gr',Rpercell',ft,opts);

para = [Ka Ke NA ksynmax NR NP Mx My rho R S k2 KMeff maa];
% Foldername = 'plotatoptimal_T3';
Foldername = 'plotatoptimal_T3_2';  %for T3 proportional to aa
plotatoptimal_T3(mfraction,odemax,para,M0f,V,M,kmetsample,Foldername,dataextract,d2,gr,Rpercell,mW_aa,ktrans_T3,Liudata,Mpercell,kselect)

%illustrate flux limited by transporter or ribosomes
%plot each kmet individually
para = [Ka Ke NA ksynmax];
for j = 1 : length(kmetsample)
%     Foldername = sprintf('kmet = %.2f_T3',kmetsample(j)/3600);
    Foldername = sprintf('kmet = %.2f_T3_2',kmetsample(j)/3600);  %for T3 proportional to aa
    eachkmet_plotting_T3(odestore,para,j,V,M,kmetsample,Foldername,ktrans_T3)
end

%check time evolution of ode for each condition
para2 = [ksynmax nSim NR NP Ke NA];
for j = 1 : length(kmetsample)
%     Foldername = sprintf('kmet = %.2f_timeode_T3',kmetsample(j)/3600);
    Foldername = sprintf('kmet = %.2f_timeode_T3_2',kmetsample(j)/3600); %for T3 proportional to aa
    recheckode_T3(j,odemax,kmetsample,flist,nstore,para2,M,Rmfc,Pmfc,aamfc,V,ktrans_T3,Foldername,Vf,T3convertf_T3)
end

for j = 1 : 3
    store = odestore{j};
    figure
    h = area(store(:,7), [store(:,7).*M(store) store(:,8).*M(store) store(:,9).*M(store)])
    xlabel('R mass fraction')
    ylabel('Cell mass (aa)')
    legend('R','P','aa')
    h(3).FaceColor = [0.4940 0.1840 0.5560];
%     title(sprintf('k_{met} = %.1f 1/s,m_{th} = %.2g aa',kmetsample(j)/3600,nstore(j,6)))
    title(sprintf('k_{met} = %.2f 1/s,m_{0} = %.2g aa',kmetsample(j)/3600,nstore(j,6)))

end

figure
for j = 1 : length(kmetsample)
    plot(odestore{j}(:,7), M(odestore{j}))
    hold on
end
scatter(mfraction(:,1),odemax(:,10),'filled')
xlabel('R mass fraction')
ylabel('Cell mass (aa)')
ylim([1e9 5e9])


figure
for j = 1 : length(kmetsample)
%     plot(uRlist, odestore{j}(:,4))
    plot(odestore{j}(:,7), odestore{j}(:,4))
    hold on
end
% scatter(odemax(:,5),odemax(:,4),'filled')
scatter(mfraction(:,1),odemax(:,4),'filled')

ylabel('Growth rate (1/h)')
xlabel('R mass fraction')

%plot all kmet in a same plot
%plot aa vs R mass fraction
figure
for j = 1 : length(kmetsample)
    semilogy(odestore{j}(:,7), odestore{j}(:,3))
    hold on
end
scatter(odemax(:,7),odemax(:,3),'filled')
xlabel('R mass fraction')
ylabel('#aa')

%fluxes and ribosome mass fraction
figure
for j = 1 : length(kmetsample)
    JMstore = odestore{j}(:,2).*kmetsample(j)./(1+(odestore{j}(:,3)./(Ke*V(odestore{j})*1e-15*NA)).^2);
    JSstore =  odestore{j}(:,1).*ksynmax.*odestore{j}(:,3)./(Ka*V(odestore{j})*NA*1e-15+odestore{j}(:,3));
    plot(odestore{j}(:,7), JMstore)
    hold on
end
scatter(odemax(:,7), odemax(:,2).*kmetsample./(1+(odemax(:,3)./(Ke*V(odemax)*1e-15*NA)).^2),'filled')    
xlabel('R mass fraction')
ylabel('Flux (aa/h)')

%to illustrate max growth rate doesn't require max flux
for j = 1 : length(kmetsample)
   JMstore = odestore{j}(:,2).*kmetsample(j)./(1+(odestore{j}(:,3)./(Ke*V(odestore{j})*1e-15*NA)).^2);
   [val, id] = max(odestore{j}(:,4));
   [val2, id2] = max(JMstore);

   figure
   yyaxis left
   plot(odestore{j}(:,7), odestore{j}(:,4))   
   hold on
   plot(odestore{j}(id,7)*ones(length(odestore{j}(:,7)),1),odestore{j}(:,4))
   ylabel('Growth rate (1/h)')
   yyaxis right
   plot(odestore{j}(:,7), JMstore)
   hold on
   plot(odestore{j}(id2,7)*ones(length(odestore{j}(:,7)),1),JMstore)
   ylabel('Intake flux (aa/h)')
   title(sprintf('k_{met} = %.2f 1/s,m_{0} = %.2g aa',kmetsample(j)/3600,nstore(j,6)))

end

%to illustrate input fluxes drive growth
% para = [Ka Ke NA ksynmax];
% fluxdrive_allkmet(odestore,para,V,M,kmetsample,odemax,nstore)

%to illustrate #R \propto gr^2 and J_M \propto M^2
% para = [Ka Ke NA ksynmax NR NP];
% plotatoptimal(mfraction,odemax,para,M0f,V,M,kmetsample)


figure
for j =  1 : length(kmetsample)
    JMstore = odestore{j}(:,2).*kmetsample(j)./(1+(odestore{j}(:,3)./(Ke*V(odestore{j})*1e-15*NA)).^2);
    [val, id] = max(odestore{j}(:,4)); %max growth rate
    frange = find(abs(JMstore-JMstore(id))<=0.3*JMstore(id));
    plot(JMstore(frange),odestore{j}(frange,4))
    hold on
end
scatter(odemax(:,2).*kmetsample./(1+(odemax(:,3)./(Ke*V(odemax)*1e-15*NA)).^2),odemax(:,4),'filled')
xlabel('Intake flux (aa/h)')
ylabel('Growth rate (1/h)')

figure
for j =  1 : length(kmetsample)
    JMstore = odestore{j}(:,2).*kmetsample(j)./(1+(odestore{j}(:,3)./(Ke*V(odestore{j})*1e-15*NA)).^2);
    JSstore =  odestore{j}(:,1).*ksynmax.*odestore{j}(:,3)./(Ka*V(odestore{j})*NA*1e-15+odestore{j}(:,3));

    [val, id] = max(odestore{j}(:,4)); %max growth rate
    frange = find(abs(JMstore-JMstore(id))<=0.3*JMstore(id));
    %     plot(JMstore, odestore{j}(:,4))
%     plot(odestore{j}(frange,7), JMstore(frange))
    plot(odestore{j}(frange,7),[JMstore(frange) JSstore(frange)])
    hold on
    plot(odestore{j}(frange,7),[odestore{j}(frange,2).*kmetsample(j) odestore{j}(frange,1).*ksynmax],'--')

end
scatter(odemax(:,7),odemax(:,2).*kmetsample./(1+(odemax(:,3)./(Ke*V(odemax)*1e-15*NA)).^2),'filled')
ylabel('Flux (aa/h)')
xlabel('R mass fraction')
ylim([0 12e9])
xlim([0 0.4])

figure
for j =  1 : length(kmetsample)
    JMstore = odestore{j}(:,2).*kmetsample(j)./(1+(odestore{j}(:,3)./(Ke*V(odestore{j})*1e-15*NA)).^2);
    [val, id] = max(odestore{j}(:,4)); %max growth rate
    frange = find(abs(JMstore-JMstore(id))<=0.3*JMstore(id));
    %     plot(JMstore, odestore{j}(:,4))
%     plot(odestore{j}(frange,7), JMstore(frange))
    plot(odestore{j}(frange,7),odestore{j}(frange,4))
    hold on
end
scatter(odemax(:,7),odemax(:,4),'filled')
ylabel('Growth rate (1/h)')
xlabel('R mass fraction')

figure
for j =  1 : length(kmetsample)
    JMstore = odestore{j}(:,2).*kmetsample(j)./(1+(odestore{j}(:,3)./(Ke*V(odestore{j})*1e-15*NA)).^2);
    [val, id] = max(odestore{j}(:,4)); %max growth rate
    frange = find(abs(JMstore-JMstore(id))<=0.3*JMstore(id));
    %     plot(JMstore, odestore{j}(:,4))
%     plot(odestore{j}(frange,7), JMstore(frange))
    plot(JMstore(frange),M(odestore{j}(frange,:)))
    hold on
end
scatter(odemax(:,2).*kmetsample./(1+(odemax(:,3)./(Ke*V(odemax)*1e-15*NA)).^2),odemax(:,10),'filled')
xlabel('Input flux (aa/h)')
ylabel('Cell mass (aa)')


%illustrate flux limited by transporter or ribosomes
%plot each kmet individually
% para = [Ka Ke NA ksynmax];
% for j = 1 : length(kmetsample)
%     eachkmet_plotting(odestore,para,j,V,M,kmetsample)
% end


figure
j = 5;
JMstore = odestore{j}(:,2).*kmetsample(j)./(1+(odestore{j}(:,3)./(Ke*V(odestore{j})*1e-15*NA)).^2);
JSstore =  odestore{j}(:,1).*ksynmax.*odestore{j}(:,3)./(Ka*V(odestore{j})*NA*1e-15+odestore{j}(:,3));
JMstore./(JSstore+odestore{j}(:,4).*odestore{j}(:,3)); %check flux balance
[val, id] = max(odestore{j}(:,4)); %max growth rate
frange = find(abs(JMstore-JMstore(id))<=0.3*JMstore(id)); %find within certain range
plot(odestore{j}(frange,7),[JMstore(frange) JSstore(frange)])
hold on
plot(odestore{j}(frange,7), [odestore{j}(frange,2).*kmetsample(j) odestore{j}(frange,1).*ksynmax],'--')
title(sprintf('k{met}^{eff} = %.2f 1/s',kmetsample(j)/3600))
xlabel('R mass fraction')
ylabel('Flux (aa/h)')
legend('J_M','J_S','J_M^{max}','J_S^{max}')
xlim([0.1 0.25])
ylim([1e9 2e9])
%pick to illustrate which partition to choose via flux
[val, id] = max(odestore{j}(:,4)); %max growth rate
flux_pick = NaN*ones(3,6);
%1st - 2 - 3 - 4 - 5 col: R mass fraction - uR - JM - JS - JMmax - JSmax
rgap = [-5 0 5];
for i = 1 : 3
   flux_pick(i,:) =  [odestore{j}(id-rgap(i),7) odestore{j}(id-rgap(i),5) JMstore(id-rgap(i)) JSstore(id-rgap(i)) odestore{j}(id-rgap(i),2).*kmetsample(j) odestore{j}(id-rgap(i),1).*ksynmax]; 
    
end

figure
bar(flux_pick(:,1),[flux_pick(:,2).*flux_pick(:,4) (1-flux_pick(:,2)).*flux_pick(:,4)],'stacked')
hold on
scatter(flux_pick(:,1),flux_pick(:,3),'filled')
plot(flux_pick(:,1),flux_pick(:,5:6),'sq--')
legend('R','P','J_M','J_M^{max}','J_S^{max}')
xlabel('R mass fraction')
ylabel('Flux (aa/h)')


%and same mass with different ribosome fraction
figure
store = odestore{j};
h = area(store(frange,7), [store(frange,7).*M(store(frange,:)) store(frange,8).*M(store(frange,:)) store(frange,9).*M(store(frange,:))])
xlabel('R mass fraction')
ylabel('Cell mass (aa)')
legend('R','P','aa')
h(3).FaceColor = [0.4940 0.1840 0.5560];
xlim([0.1 0.25])
%and growth rate with different ribosome fraction
figure
plot(odestore{j}(frange,7),odestore{j}(frange,4))
xlabel('R mass fraction')
ylabel('Growth rate (1/h)')
xlim([0.1 0.25])
%synthesis allocation vs ribosome mass fraction
figure
area(odestore{j}(frange,7),[odestore{j}(frange,5) odestore{j}(frange,6)])
xlabel('R mass fraction')
ylabel('Synthesis allocation')
legend('R','P')
xlim([0.1 0.25])
%#R,P,aa vs ribosome mass fraction
semilogy(odestore{j}(frange,7),odestore{j}(frange,1:3))
namelist = {'R','P','aa'};
for i = 1 : 3
figure
plot(odestore{j}(frange,7),odestore{j}(frange,i))
xlabel('R mass fraction')
ylabel(sprintf('#%s (molecules)',namelist{i}))
% xlim([0.1 0.25])
end

JM = odemax(:,2).*kmetsample./(1+(odemax(:,3)./(Ke*V(odemax)*1e-15*NA)).^2);
JS = odemax(:,1).*ksynmax.*odemax(:,3)./(Ka*V(odemax)*NA*1e-15+odemax(:,3));
JM./(JS+odemax(:,4).*odemax(:,3))
JM./(odemax(:,4).*odemax(:,11))
%to illustrate #R \propto gr^2 and J_M \propto M^2
% para = [Ka Ke NA ksynmax NR NP];
% plotatoptimal(mfraction,odemax,para,M0f,V,M,kmetsample)

figure
scatter(odemax(:,4),odemax(:,10),'filled')
hold on
plot(0:0.1:2,M0f(0:0.1:2))
xlabel('Growth rate (1/h)')
ylabel('Cell mass (aa)')

hfit = fit(odemax(:,4),mfraction(:,1),'poly1')
plot(hfit,odemax(:,4),mfraction(:,1))

fR = @(x) 1/NR*M0f(x).*(hfit.p1*x + hfit.p2);

figure
scatter(odemax(:,4),odemax(:,1),'filled')
hold on
plot(0:0.1:2,fR(0:0.1:2))
xlabel('Growth rate (1/h)')
ylabel('#R/cell')
set(gca,'XScale','log','YScale','log')

fit(odemax(:,4),odemax(:,1),'poly2')
%fit ksyn* vs gr* (for old values)
x0 = [ksynmax 1 7*3600];
x0 = [15 0.25 3];
fitfun = fittype(@(a,b,c,x) a*x./(b+x)+c );
[fitted_curve,gof] = fit(odemax(:,4),ksynmax/3600.*odemax(:,3)./(Ka*V(odemax)*NA*1e-15+odemax(:,3)),fitfun,'StartPoint',x0);
coeffvals = coeffvalues(fitted_curve);

figure
scatter(odemax(:,4),1/3600*ksynmax.*odemax(:,3)./(Ka*V(odemax)*NA*1e-15+odemax(:,3)),'filled')
hold on
plot(0:0.1:2,1/3600*fitted_curve(0:0.1:2))
xlabel('Growth rate (1/h)')
ylabel('Synthesis speed (1/s)')
ylim([0 22])

% fit ksyn* vs gr* for new values
ft = fittype( 'a*x/(b+x)+c', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0 0.22 0];
opts.Robust = 'Bisquare';
opts.StartPoint = [15 0.3 4];

[fitresult, gof] = fit(odemax(:,4),1/3600*ksynmax.*odemax(:,3)./(Ka*V(odemax)*NA*1e-15+odemax(:,3)), ft, opts );
ffit = @(x) 3600*(fitresult.a*x./(fitresult.b+x)+fitresult.c);

figure
scatter(odemax(:,4),1/3600*ksynmax.*odemax(:,3)./(Ka*V(odemax)*NA*1e-15+odemax(:,3)),'filled');
hold on
plot(0:0.01:2,ffit(0:0.01:2)/3600)
xlim([0 2])
ylim([0 18])
xlabel('Growth rate (1/h)')
ylabel('Synthesis speed (1/s)')


% fJMapprox = @(x) fR(x).*(coeffvals(1)*x./(coeffvals(2)+x)+coeffvals(3));
fJMapprox = @(x) fR(x).*ffit(x);
JM = odemax(:,2).*kmetsample./(1+(odemax(:,3)./(Ke*V(odemax)*1e-15*NA)).^2);

figure
scatter(odemax(:,4),JM,'filled')
hold on
plot(0:0.1:2,fJMapprox(0:0.1:2))
xlabel('Growth rate (1/h)')
ylabel('Input flux (aa/h)')
legend('ODEs','Analytical fitting')

figure
scatter(M(odemax),JM,'filled')
hold on
plot(M0f(0.2:0.1:2.5),fJMapprox(0.2:0.1:2.5))
xlabel('Cell mass (aa)')
ylabel('Input flux (aa/h)')
legend('ODEs','Analytical fitting')
set(gca,'XScale','log','YScale','log')
xlim([1e8 1e10])
ylim([1e8 1e11])

gfit = fit(log(M(odemax)),log(JM),'poly1')
% plot(gfit,log(M(maxlist)),log(JM))
figure
scatter(M(odemax),JM,'filled')
hold on
Mlist = logspace(8.5,9.5,100);
plot(Mlist, exp(gfit.p2)*Mlist.^(gfit.p1))
set(gca,'XScale','log','YScale','log')
xlabel('Cell mass (aa)')
ylabel('Input flux (aa/h)')

figure
scatter(odemax(:,4), JM,'filled')
xlabel('Growth rate (1/h)')
ylabel('Input flux (aa/h)')

figure
bar(odemax(:,4), [odemax(:,5) odemax(:,6)],'stacked')
ylabel('Synthesis allocation')
xlabel('Growth rate (1/h)')
legend('R','P')

figure
h = bar(odemax(:,4), odemax(:,7:9),'stacked')
ylabel('Mass fraction')
xlabel('Growth rate (1/h)')
legend('R','P','aa')
h(3).FaceColor = [0.4940 0.1840 0.5560];
h(2).FaceColor = [1 1 0];

figure
bar(odemax(:,4), [odemax(:,5).*JS odemax(:,6).*JS],'stacked')
hold on
scatter(odemax(:,4), JM,'filled')
ylabel('Synthesis flux (aa/h)')
xlabel('Growth rate (1/h)')
legend('R','P','J_M')

namelist = {'R','P','aa'};
masslist = [NR NP 1];

for i = 1 : length(namelist)
    yL = sprintf('#%s/cell',namelist{i}); 
    yR = sprintf('%s mass fraction',namelist{i}); 
    figure
    yyaxis left
    scatter(odemax(:,4),odemax(:,i),'filled')
    ylabel(yL)
    yyaxis right
    scatter(odemax(:,4),mfraction(:,i),'filled')
    xlabel('Growth rate (1/h)')
    ylabel(yR) 
end

figure
scatter(odemax(:,4),ksynmax*odemax(:,3)./(Ka*V(odemax)*NA*1e-15+odemax(:,3))*1/3600,'filled')
xlabel('Growth rate (1/h)')
ylabel('Synthesis speed (1/s)')

figure
scatter(odemax(:,4),kmetsample*1/3600,'filled')
xlabel('Growth rate (1/h)')
ylabel('k^{eff}_{met} (1/s)')

figure
scatter(kmetsample.*0.05./(0.4*3600-kmetsample),odemax(:,4),'filled')
ylabel('Growth rate (1/h)')
xlabel('[Nutrients] (mM)')
%fit kmet* vs gr*
x0 = [0.2*3600 0.01];
fitfun = fittype(@(a,b,x) a*x./(b+x) );
[fitted_curve,gof] = fit(kmetsample.*0.05./(0.4*3600-kmetsample),kmetsample,fitfun,'StartPoint',x0);
coeffvals_kmet = coeffvalues(fitted_curve);
figure
scatter(kmetsample.*0.05./(0.4*3600-kmetsample),kmetsample/3600,'filled')
hold on
plot(0:0.01:0.2,fitted_curve(0:0.01:0.2)/3600)
ylabel('k^{eff}_{met} (1/s)')
xlabel('[Nutrients] (mM)')

figure
scatter(odemax(:,4), M(odemax),'filled')
xlabel('Growth rate (1/h)')
ylabel('Cell mass (aa)')

%% one kmet, scan ribosome synthesis allocation uR to select fastest growth rate
j = 3;
kmet = kmetsample(j);
store = odestore{j};
[val, id] = max(store(:,4));

figure
plot(store(:,7), store(:,4))
hold on
scatter(store(id,7), store(id,4),'filled')
xlabel('R mass fraction')
ylabel('Growth rate (1/h)')

namelist = {'R','P','aa'};
masslist = [NR NP 1];
for i = 1 : length(namelist)
    figure
    yyaxis left
    semilogy(store(:,4), store(:,i))
    hold on
    scatter(store(id,4), store(id,i),'filled')
    ylabel(sprintf('#%s',namelist{i}))

    yyaxis right
    plot(store(:,4), masslist(i)*1e6*store(:,i)./(V(store)*1e-15*NA))
    hold on
    scatter(store(id,4), masslist(i)*1e6*store(id,i)./(Vf(store(id,:))*1e-15*NA),'filled')
    set(gca,'Yscale','log')
    ylabel(sprintf('[%s] (\\muM)',namelist{i}))
%     xlabel('R mass fraction')  
    xlabel('Growth rate (1/h)')

end

for i = 1 : length(namelist)
    figure
    yyaxis left
    semilogy(store(:,4), store(:,i))
    hold on
    scatter(store(id,4), store(id,i),'filled')
    ylabel(sprintf('#%s',namelist{i}))

    yyaxis right
    plot(store(:,4), masslist(i)*store(:,i)./(M(store)))
    hold on
    scatter(store(id,4), masslist(i)*store(id,i)./(Mf(store(id,:))),'filled')
    ylabel(sprintf('%s mass fraction',namelist{i}))
%     xlabel('R mass fraction')
    xlabel('Growth rate (1/h)')

end

%input flux 
JM = store(:,2).*kmet./(1+(store(:,3)./(Ke*V(store)*NA*1e-15)).^2);
% JM_end = store_end(:,2).*kmet./(1+(store_end(:,3)./(Ke*V(store_end)*NA*1e-15)).^2);

%synthesis flux 
ktrans = ktrans_T3;
JS = store(:,1).*ktrans(store);
% JS_end = store_end(:,1).*ktrans(store_end);


% figure
% plot(store(:,7), fluxave(:,1)) %fluxave = JM/log(2) 
% hold on
% plot(store(:,7), [JM JS])
% xlabel('R mass fraction')
% ylabel('Flux (aa/h)')
% legend('J_M','J_S','J^{max}_M','J^{max}_S')

figure
plot(store(:,7), [JM JS])
hold on
plot(store(:,7), [store(:,2)*kmet store(:,1)*ksynmax],'--')
xlabel('R mass fraction')
ylabel('Flux (aa/h)')
legend('J_M','J_S','J^{max}_M','J^{max}_S')

figure
area(store(:,7), [store(:,5).*JS store(:,6).*JS])
xlabel('R mass fraction')
ylabel('Synthesis flux (aa/h)')
legend('R','P')
xlim([0 0.4])

figure
h = area(store(:,7), [store(:,7).*M(store) store(:,8).*M(store) store(:,9).*M(store)])
xlabel('R mass fraction')
ylabel('Cell mass (aa)')
legend('R','P','aa')
h(3).FaceColor = [0.4940 0.1840 0.5560];

figure
h = area(store(:,7), [store(:,7) store(:,8) store(:,9)])
xlabel('R mass fraction')
ylabel('Mass fraction')
legend('R','P','aa')
h(3).FaceColor = [0.4940 0.1840 0.5560];

figure
area(store(:,7), [store(:,5) store(:,6)])
xlabel('R mass fraction')
ylabel('Synthesis allocation')
legend('R','P')
