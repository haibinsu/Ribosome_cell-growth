load('Lidata')
%from https://doi.org/10.1038/s41564-020-0717-x
%1st - 2nd : growth rate (1/h) - mass per cell (OD600 ml per 1e9 cells)
id = find(Lidata(:,1) >= 0);
[mfit,gof] = fit(Lidata(id,1),Lidata(id,2),'poly1');
figure
scatter(Lidata(:,1),Lidata(:,2),'filled')
hold on
plot(0.05:0.01:1.8,mfit(0.05:0.01:1.8))
xlabel('Growth rate (1/h)')
ylabel('Average cell mass (OD600 ml per 10^9 cells)')
saveas(gca,'Mass_gr.eps','eps')


load('growthdata')
load('growthCm')
%1st - 4th = growth rate, RNA mass/OD (ug), total protein mass/OD (ug), #cell/OD600 
%5th col: cell volume (um^3)
[~,id] = sort(growthdata(:,1)); %sort growth rate from slow to fast 
growthdata = growthdata(id,:);
[~,id] = sort(growthCm(:,1)); %sort growth rate from slow to fast 
growthCm = growthCm(id,:);

%compare protein density
figure
plot(growthCm(:,1),growthCm(:,3)./(growthCm(:,4).*growthCm(:,5))/1e-8,'.-')
hold on
plot(growthdata(:,1),growthdata(:,3)./(growthdata(:,4).*growthdata(:,5))/1e-8,'sq-')
ylim([8 18])
ylabel('Total protein density (10^{-8} \mug/\mum^3)')
xlabel('Growth rate (1/h)')
legend('chloramphenicol','nutrient controlled')
saveas(gca,'density_gr.eps','eps')

%in Bremmer and Dennis  10.1128/ecosal.5.2.3
gr = [0.6 1 1.5 2 2.5 3]*log(2); %growth rate 1/h
RNApercell = [23 44 76 128 180 214]; %ug/10^9 cells
Rpercell = [8 15 26 44 61 73]*1e3; 
x0 = [1 1];
fitfun = fittype( @(a,b,x) a*x.^2+b*x);
[fitted_curve,gof] = fit(gr',RNApercell',fitfun,'StartPoint',x0);
coeffvals = coeffvalues(fitted_curve);
[hfit,gof] = fit(gr',RNApercell','exp1');
[gfit,gof] = fit(gr',Rpercell',fitfun,'StartPoint',x0);

figure
scatter(gr,RNApercell, 'filled')
hold on
plot(gr,fitted_curve(gr))
plot(gr,hfit(gr))
hold off
xlabel('Growth rate (1/h)')
ylabel('RNA/cell (\mug/10^9 cells)')
% ylabel('R/cell')
legend('Experiments','Fit ax^2+bx','Fit ae^{bx}')

figure
scatter(gr,Rpercell, 'filled')
hold on
plot(0:0.1:2,gfit(0:0.1:2))
hold off
xlabel('Growth rate (1/h)')
ylabel('R/cell')
legend('Experiments','Fit ax^2+bx')
xlim([0 2.5])

%compare protein density
% figure
% yyaxis left
% plot(growthCm(:,1),growthCm(:,3)./(growthCm(:,4).*growthCm(:,5))/1e-8,'.-')
% hold on
% plot(growthdata(:,1),growthdata(:,3)./(growthdata(:,4).*growthdata(:,5))/1e-8,'sq-')
% ylim([0 20])
% ylabel('Total protein density (10^{-8} \mug/\mum^3)')
% yyaxis right
% Lyax = gca;  
% Lyaxt = Lyax.YAxis(1).TickValues
% Lyaxconvert = Lyaxt*1e-8*1e-6/(110*1.66e-24);     % Convert to #aa 1 aa = 110*1.66e-24 g
% Lyax.YAxis(2).TickLabelFormat = '%e'
% Lyax.YAxis(2).TickLabels = Lyaxconvert;
% ylabel('Total protein density (aa/\mum^3)')
% xlabel('Growth rate (1/h)')
% legend('chloramphenicol','nutrient controlled')


%from Hwa
fitfun = fittype( @(a,b,x) a*x.^2+b*x );
[gfit,gof] = fit(growthdata(:,1),(growthdata(:,2)./growthdata(:,4))/1e-9,fitfun,'StartPoint',x0);
coeffvals = coeffvalues(gfit);
[mfit,gof] = fit(growthdata(:,1),(growthdata(:,2)./growthdata(:,4))/1e-9,'exp1');

figure
scatter(growthdata(:,1),(growthdata(:,2)./growthdata(:,4))/1e-9, 'filled')
hold on
plot(growthdata(:,1),gfit(growthdata(:,1)))
plot(growthdata(:,1),mfit(growthdata(:,1)))
hold off
xlabel('Growth rate (1/h)')
ylabel('RNA/cell (\mug/10^9 cells)')
legend('Experiments','Fit ax^2+bx','Fit ae^{bx}')


%glucose intake flux per 1 gDW vs growth rate
scatter([0.4 0.7 1 1.27 1.84], [4 8 10.5 12 12.2])
xlabel('Growth rate (1/h)')
ylabel('Intake flux (mmol/gDW/h)')
%average cell dry mass vs growth rate
scatter([0.4 0.7 1 1.27 1.84], [0.27 0.33 0.47 0.63 1.48])   
xlabel('Growth rate (1/h)')
ylabel('Average cell DW (pg)')
%%glucose intake flux per cell vs growth rate
scatter([0.4 0.7 1 1.27 1.84], 1e-3*6.02e23*[4 8 10.5 12 12.2].*([0.27 0.33 0.47 0.63 1.48]*1e-12) )
xlabel('Growth rate (1/h)')
ylabel('Intake flux (#glucose/cell/h)')

%nutrient 
%volume per cell vs gr
scatter(growthdata(:,1),growthdata(:,5))
xlabel('Growth rate (1/h)')
ylabel('Cell volume (\mum^3) ')

%derivative
plot(growthdata(1:end-1,5),diff(growthdata(:,5))./diff(growthdata(:,1)),'.-')
ylabel('dV/demp\lambda (\mug.h)')
xlabel('Cell volume (\mum^3) ')

%RNA mass fraction 
scatter(growthdata(:,1),growthdata(:,2)./growthdata(:,3))
xlabel('Growth rate (1/h)')
ylabel('\phi_R Ribosome mass fraction ')

%derivative
temp2 = growthdata(:,2)./growthdata(:,3);
plot(temp2(1:end-1),diff(temp2)./diff(growthdata(:,1)),'.-')
ylabel('d\phi_R/d\lambda (h)')
xlabel('\phi_R ')

%dV/dgr vs dphiR/dgr
plot(diff(growthdata(:,5))./diff(growthdata(:,1)),diff(temp2)./diff(growthdata(:,1)),'.-')
ylabel('d\phi_R/d\lambda (h)')
xlabel('dm/d\lambda (\mug.h)')


%antibiotics 
%volume per cell
scatter(growthCm(:,1),growthCm(:,5))
xlabel('Growth rate (1/h)')
ylabel('Cell volume (\mum^3) ')
title('chloramphenicol')

%derivative
plot(growthCm(1:end-1,5),diff(growthCm(:,5))./diff(growthCm(:,1)),'.-')
ylabel('dV/d\lambda (\mum^3.h)')
xlabel('Cell volume (\mum^3) ')
title('chloramphenicol')

%RNA mass fraction 
scatter(growthCm(:,1),growthCm(:,2)./growthCm(:,3))
xlabel('Growth rate (1/h)')
ylabel('\phi_R Ribosome mass fraction ')
title('chloramphenicol')

%derivative
temp2 = growthCm(:,2)./growthCm(:,3);
plot(temp2(1:end-1),diff(temp2)./diff(growthCm(:,1)),'.-')
ylabel('d\phi_R/d\lambda (h)')
xlabel('\phi_R ')
title('chloramphenicol')

%dV/dgr vs dphiR/dgr
plot(diff(growthCm(:,5))./diff(growthCm(:,1)),diff(temp2)./diff(growthCm(:,1)),'.-')
ylabel('d\phi_R/d\lambda (h)')
xlabel('dV/d\lambda (\mum^3.h)')
title('chloramphenicol')

%protein mass/cell and cell volume (Cm)
yyaxis left
plot(growthCm(:,1),(growthCm(:,3)./growthCm(:,4))/1e-8,'.-')
ylabel('mass per cell(10^{-8}\mug)')

yyaxis right
plot(growthCm(:,1),growthCm(:,5),'o-')
ylabel('cell volume (\mum^3)')
xlabel('growth rate (1/h)')

%RNA/cell

plot(growthCm(:,1),(growthCm(:,2)./growthCm(:,4))/1e-8,'.-')
ylabel('RNA mass per cell(10^{-8}\mug)')
xlabel('growth rate (1/h)')

% plot(growthCm(:,1),(growthCm(:,3)./growthCm(:,4)),'.-')
%other protein/cell
plot(growthCm(:,1),(growthCm(:,3)-growthCm(:,2))./growthCm(:,4)*1/1e-8,'.-')
ylabel('other protein mass per cell(10^{-8}\mug)')
xlabel('growth rate (1/h)')


%compare protein density
figure
plot(growthCm(:,1),growthCm(:,3)./(growthCm(:,4).*growthCm(:,5))/1e-8,'.-')
hold on
plot(growthdata(:,1),growthdata(:,3)./(growthdata(:,4).*growthdata(:,5))/1e-8,'sq-')
ylim([8 18])
ylabel('Total protein density (10^{-8} \mug/\mum^3)')
xlabel('Growth rate (1/h)')
legend('chloramphenicol','nutrient controlled')

%protein mass/cell and cell volume (nutrient controlled)
yyaxis left
plot(growthdata(:,1),(growthdata(:,3)./growthdata(:,4))/1e-8,'.-')
ylabel('mass per cell(10^{-8}\mug)')

yyaxis right
plot(growthdata(:,1),growthdata(:,5),'o-')
ylabel('cell volume (\mum^3)')
xlabel('growth rate (1/h)')

%RNA/cell

plot(growthdata(:,1),(growthdata(:,2)./growthdata(:,4))/1e-8,'.-')
ylabel('RNA mass per cell(10^{-8}\mug)')
xlabel('growth rate (1/h)')

% plot(growthCm(:,1),(growthCm(:,3)./growthCm(:,4)),'.-')
%other protein/cell
plot(growthdata(:,1),(growthdata(:,3)-growthdata(:,2))./growthdata(:,4)*1/1e-8,'.-')
ylabel('other protein mass per cell(10^{-8}\mug)')
xlabel('growth rate (1/h)')

hfit = fit(growthdata(:,1),growthdata(:,3)./growthdata(:,4),'exp1')
M0 = hfit.a;
aM = hfit.b;

plot(hfit,growthdata(:,1),growthdata(:,3)./growthdata(:,4))
xlabel('Growth rate (1/h)')
ylabel('Mass per cell (\mug)')
%RNA mass per cell 
gfit = fit(growthdata(:,1),growthdata(:,2)./growthdata(:,4),'exp1')
R0 = gfit.a;
bR = gfit.b;

plot(gfit,growthdata(:,1),growthdata(:,2)./growthdata(:,4))
xlabel('Growth rate (1/h)')
ylabel('RNA per cell (\mug)')
%mass- growth rate
g = @(x) M0*exp(aM*x);

%RNA - growth rate
h = @(x) R0*exp(bR*x);

%mass fraction of ribosome 
f = @(x) R0/M0*exp((bR-aM)*x);

%linear fit RNA mass fraction
ufit = fit(growthdata(:,1),growthdata(:,2)./growthdata(:,3),'poly1')
p1R = ufit.p1;
p2R = ufit.p2;


%RNA mass fraction 
scatter(growthdata(:,1),growthdata(:,2)./growthdata(:,3))
hold on
plot(gr,f(gr))
plot(gr,p1R*gr+p2R)
xlabel('Growth rate (1/h)')
ylabel('\phi_R Ribosome mass fraction ')

legend('Exp','exponential fit ratio','linear fit')

scatter(growthdata(:,1),(growthdata(:,3)-growthdata(:,2))./growthdata(:,4))
hold on
plot(gr,g(gr) - h(gr))
xlabel('Growth rate (1/h)')
ylabel('Remaining protein mass (\mug) ')

NP= 400;
Ymax = 0.45; %gDW/gC
Ymax*1/maa/(1/180.16*NA)

scatter(growthdata(:,1),1/(Ymax*1/maa/(1/180.16*NA))*NP*growthdata(:,1)./(1-growthdata(:,2)./growthdata(:,3)))
hold on
plot(gr,1/(Ymax*1/maa/(1/180.16*NA))*NP*gr./(1-h(gr)./g(gr)))
%run uptake_glucose first 
errorbar(store_clean(:,1),store_clean(:,2)*1e-3*NA/(1/(maa*NP)),qpos*1e-3*NA/(1/(NP*maa)),qneg*1e-3*NA/(1/(NP*maa)),grpos,grneg)

xlabel('Growth rate (1/h)')
ylabel('Intake rate (1/h)')
legend('Exp','Fit','exp uptake data')

scatter(mW(:,1),3600*ratio*mW(:,2)./(mW(:,1)/(NP*maa)))
hold on
scatter(growthdata(:,3)./growthdata(:,4)*1e-6,1/(Ymax*1/maa/(1/180.16*NA))*NP*growthdata(:,1)./(1-growthdata(:,2)./growthdata(:,3)),'filled')
set(gca,'YScale','log','XScale','log')
xlabel('Cell mass (g)')
ylabel('Intake rate (#aa/(#P.h)')



phiR = 0:0.01:0.5;

plot(phiR,phiR.*22.*phiR./(0.11+phiR))


scatter(growthCm(:,1),growthCm(:,2)./growthCm(:,4))
scatter(growthCm(:,1),growthCm(:,3)./growthCm(:,4))
scatter(growthCm(:,1),growthCm(:,2)./growthCm(:,3))

mfit = fit(growthCm(:,1),growthCm(:,3)./growthCm(:,4),'exp1')
plot(mfit,growthCm(:,1),growthCm(:,3)./growthCm(:,4))