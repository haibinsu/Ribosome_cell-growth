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


%use this one 
Mpick = [1.14075e9 1.6476e9 1.84e9 2.2691334e9 2.6478e9];
kmetsample = [0.07 0.1 0.14 0.17 0.24]'*3600;
nstore = NaN*ones(length(kmetsample),6);
%1 - 6th col: #R - #P - #aa - gr1 - gr2 - M0

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

%to get R mass fraction vs growth rate 
%to illustrate input fluxes drive growth
para = [Ka Ke NA ksynmax];
% Foldername = 'allkmetsameplot_T3';
Foldername = 'allkmetsameplot_T3_2'; %for T3 proportional to aa
fluxdrive_allkmetT3(odestore,para,V,M,kmetsample,odemax,nstore,Foldername,ktrans_T3,glawdata)

%to get #R/cell vs growth rate and JM vs growth rate 
%and gr/JM vs gr 
%in Bremmer and Dennis  10.1128/ecosal.5.2.3
gr = [0.6 1 1.5 2 2.5 3]*log(2); %growth rate 1/h
RNApercell = [23 44 76 128 180 214]; %ug/10^9 cells
Rpercell = [8 15 26 44 61 73]*1e3; 
kselect = [13 18 21 22 22 22]*3600; %aa/(R.h)

para = [Ka Ke NA ksynmax NR NP Mx My rho R S k2 KMeff maa];
Foldername = 'plotatoptimal_T3_2';  %for T3 proportional to aa
plotatoptimal_T3(mfraction,odemax,para,M0f,V,M,kmetsample,Foldername,dataextract,d2,gr,Rpercell,mW_aa,ktrans_T3,Liudata,Mpercell,kselect)
