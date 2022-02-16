function plotatoptimal_T3(mfraction,odemax,para,M0f,V,M,kmetsample,Foldername,dataextract,d2,gr,Rpercell,mW_aa,ktrans_T3,Liudata,Mpercell,kselect)
%use KMeff from part 1 for ternary complex instead of rescaling 

Ka = para(1);
Ke = para(2);
NA = para(3);
ksynmax = para(4);
NR = para(5);
NP = para(6);
Mx = para(7);
My = para(8);
rho = para(9);
R = para(10);
S = para(11);
k2 = para(12);
KMeff = para(13);
maa = para(14);
%illustrate flux limited by transporter or ribosomes
% fpath = 'D:\Onedrive-NTU\OneDrive - Nanyang Technological University\PhD\Kinetics\Kinetics methods\OPTI-master\Figures\growth_HB\';
% yourFolder = [fpath,sprintf('kmet = %.2f',kmet/3600)];
% yourFolder = ['../OPTI-master/Figures/growth_HB/' 'plotatoptimal_T3'];
yourFolder = ['../growth_HB/' Foldername];

if exist(yourFolder, 'dir') ~= 7 %folder does not exist
       mkdir(yourFolder)
end

%to illustrate #R \propto gr^2 and J_M \propto M^2
% c = [0 0 0; 96 60 38; 142 88 56; 193 151 96; 255 199 127]/255; %list of color
c = [77 191 237; 113 152 241; 148 114 244; 188 71 248; 255 0 255]/255; %list of color
figure
plot(0.4:0.1:2,M0f(0.4:0.1:2),'LineWidth',1.2,'Color','k')
hold on
scatter(odemax(:,4),odemax(:,11),250,c,'filled')
scatter(gr,Mpercell*0.6,290,'kx','LineWidth',1.5)  
% RNA+Protein = 0.7cell dry mass, RNA = M_R/0.76 SI of 10.1126/science.1192588
% --> mass ribosome + protein ~ 0.5 - 0.6 cell dry mass
% scatter(gr,Mpercell,50,'k','sq')  
% scatter(Liudata(:,1),Liudata(:,2)*1e-12/maa,'filled')
xlabel('\lambda^* (1/h)')
ylabel('M^* (aa)')
set(gca,'FontSize',30)
xticks(0:0.4:2)
xlim([0 2.2])
box on
set(gca,'FontName','Helvetica')
saveas(gca,fullfile(yourFolder,'mass_gr'),'png')
%for matlab at home 
%solution based on https://stackoverflow.com/questions/10985946/how-to-export-the-figure-to-color-eps-in-matlab
saveas(gca,fullfile(yourFolder,'mass_gr.eps'),'epsc')

hfit = fit(odemax(:,4),mfraction(:,1),'poly1')
plot(hfit,odemax(:,4),mfraction(:,1))

%gr vs #R/cell 
fR = @(x) 1/NR*M0f(x).*(hfit.p1*x + hfit.p2);
%Hwa data fitting
ft = fittype( 'a*x^2+b*x+c', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [1e4 1000 100];
opts.Robust = 'Bisquare';
opts.StartPoint = [1e4 5000 300];

[fitresult, gof] = fit(dataextract(:,3),dataextract(:,1), ft, opts );
% ffit = @(x) 3600*(fitresult.a*x./(fitresult.b+x)+fitresult.c);
ffit = @(x) fitresult.a*x.^2+fitresult.b*x+fitresult.c;

%Bremer-Dennis fitting
ft = fittype( 'a*x^2+b*x+c', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [1e3 100 1000];
opts.Robust = 'Bisquare';
opts.StartPoint = [1e4 5000 300];
[gfit,gof] = fit(gr',Rpercell',ft,opts);

figure
plot(0:0.1:2,fR(0:0.1:2),'k')
hold on
scatter(odemax(:,4),odemax(:,1),50,c,'filled')
scatter(dataextract(:,3),dataextract(:,1),50,'sq','MarkerEdgeColor','k','LineWidth',1.25,'MarkerEdgeAlpha',1) %[0 0.5 0] [0.49 0.18 0.56]
scatter(gr,Rpercell,50,'x','MarkerEdgeColor','k','LineWidth',1.25,'MarkerEdgeAlpha',1) %[0.87 0.49 0]
xlabel('Growth rate (1/h)')
ylabel('#R/cell')
legend('Theory','Numerical results','Dai et al.','Bremer & Dennis','Location','NorthWest')
ylim([0 10e4])
yticks([0:1e4:10e4]) 
xticks(0:0.2:2) 
xlim([0 2.1])
set(gca,'FontName','Helvetica','FontSize',16)
saveas(gca,fullfile(yourFolder,'R_gr'),'png')
saveas(gca,fullfile(yourFolder,'R_gr.eps'),'epsc')
saveas(gca,fullfile(yourFolder,'R_gr.pdf'),'pdf')

figure
plot(fR(0:0.1:2),0:0.1:2,'k')
hold on
scatter(odemax(:,1),odemax(:,4),50,c,'filled')
scatter(dataextract(:,1),dataextract(:,3),50,'sq','MarkerEdgeColor','k','LineWidth',1.25,'MarkerEdgeAlpha',1) %[0 0.5 0] [0.49 0.18 0.56]
scatter(Rpercell,gr,50,'x','MarkerEdgeColor','k','LineWidth',1.25,'MarkerEdgeAlpha',1) %[0.87 0.49 0]
ylabel('\lambda^* (1/h)')
xlabel('R^* (molecules/cell)')
% legend('Theory','Numerical results','Dai et al.','Bremer & Dennis','Location','NorthWest')
xlim([0 10e4])
xticks([0:2e4:10e4]) 
yticks(0:0.4:2) 
ylim([0 2.1])
get(gca,'fontname')
set(gca,'FontName','Helvetica','FontSize',16)
saveas(gca,fullfile(yourFolder,'R_gr_flip'),'png')
saveas(gca,fullfile(yourFolder,'R_gr_flip.eps'),'epsc')
saveas(gca,fullfile(yourFolder,'R_gr_flip.pdf'),'pdf')


figure
scatter(odemax(:,4),odemax(:,10),50,c,'filled')
xticks(0:0.2:2) 
xlim([0 2])
xlabel('Growth rate (1/h)')
ylabel('#T_3/cell')
set(gca,'FontName','Helvetica')
saveas(gca,fullfile(yourFolder,'T3_gr'),'png')
saveas(gca,fullfile(yourFolder,'T3_gr.eps'),'epsc')

% figure
% plot(0:0.1:2,fR(0:0.1:2))
% hold on
% scatter(odemax(:,4),odemax(:,1),[],c,'filled')
% plot(0:0.01:2,ffit(0:0.01:2),'Color',[0.49 0.18 0.56])
% scatter(dataextract(:,3),dataextract(:,1),'sq','filled','MarkerFaceColor',[0.49 0.18 0.56])
% plot(0:0.01:2.1,gfit(0:0.01:2.1),'Color',[0.47 0.67 0.19])
% scatter(gr,Rpercell,'diamond','filled','MarkerFaceColor',[0.47 0.67 0.19])
% xlabel('Growth rate (1/h)')
% ylabel('#R/cell')
% legend('Theory','Modelling','Fitting','Dai et al.','Fitting','Bremer & Dennis','Location','NorthWest')
% ylim([0 14e4])
% xticks(0:0.2:2) 
% xlim([0 2.1])
% % set(gca,'XScale','log','YScale','log')
% saveas(gca,fullfile(yourFolder,'R_gr'),'png')
% saveas(gca,fullfile(yourFolder,'R_gr.eps'),'epsc')

ft = fittype( '(a*x+c)/(b+x)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [10 0.22 1];
opts.Upper = [22 1 5];
opts.Robust = 'Bisquare';
opts.StartPoint = [10 0.3 4];

[fitresult, gof] = fit(odemax(:,4),1/3600*ktrans_T3(odemax), ft, opts );
ffit = @(x) 3600*(fitresult.a*x+fitresult.c)./(fitresult.b+x);
ftheoryapprox = @(x) (NR/(hfit.p1)*x+NR*My/(Mx*hfit.p1))./((hfit.p1*My+hfit.p2*Mx)/(Mx*hfit.p1)+x);
ftheory = @(x) x.*M0f(x)./fR(x);

figure
% plot(0:0.01:2,ffit(0:0.01:2)/3600,'k')
plot(0:0.01:2,ftheory(0:0.01:2)/3600,'k')
hold on
plot(0.4:0.01:2,ftheoryapprox(0.4:0.01:2)/3600,'k--')
scatter(odemax(:,4),1/3600*ktrans_T3(odemax),50,c,'filled');
% scatter(data_Dai(:,1),data_Dai(:,2),50,'filled');
legend('Theory','Approximation','Location','SouthEast')
xlim([0 2])
xticks(0:0.2:2)
ylim([0 20])
xlabel('Growth rate (1/h)')
ylabel('Synthesis speed (1/s)')
saveas(gca,fullfile(yourFolder,'ksyn_gr'),'png')
saveas(gca,fullfile(yourFolder,'ksyn_gr.eps'),'epsc')


% fJMapprox = @(x) fR(x).*(coeffvals(1)*x./(coeffvals(2)+x)+coeffvals(3));
fJMapprox = @(x) fR(x).*ffit(x);
% fJMapprox = @(x) fR(x).*ftheory(x) + 0.03*M0f(x); %x is growth rate
JM = odemax(:,2).*kmetsample./(1+(odemax(:,3)./(Ke*V(odemax)*1e-15*NA)).^2);
JS = odemax(:,1).*ktrans_T3(odemax);

figure
% c = [0 0 0; 96 60 38; 142 88 56; 193 151 96; 255 199 127]/255; %list of color
c = [77 191 237; 113 152 241; 148 114 244; 188 71 248; 255 0 255]/255; %list of color
scatter(Rpercell.*kselect,gr,50,'kx')
hold on
scatter(dataextract(:,1).*dataextract(:,4)*3600,dataextract(:,3),50,'ko')
plot(fJMapprox(0:0.1:2),0:0.1:2,'k')
hold on
scatter(JM,odemax(:,4),50,c,'filled')
ylabel('Growth rate (1/h)')
xlabel('Intake flux (aa/h)')
box on
set(gca,'FontName','Helvetica')
legend('Bremer-Dennis','Hwa','Theory line','Numerical results')
saveas(gca,fullfile(yourFolder,'Intakeflux_gr'),'png')
%for matlab at home 
%solution based on https://stackoverflow.com/questions/10985946/how-to-export-the-figure-to-color-eps-in-matlab
saveas(gca,fullfile(yourFolder,'Intakeflux_gr.eps'),'epsc')

figure
% c = [0 0 0; 96 60 38; 142 88 56; 193 151 96; 255 199 127]/255; %list of color
c = [77 191 237; 113 152 241; 148 114 244; 188 71 248; 255 0 255]/255; %list of color
yyaxis left
plot(fJMapprox(0:0.1:2),0:0.1:2)
hold on
scatter(JM,odemax(:,4),50,c,'filled')
ylabel('Growth rate (1/h)')
yyaxis right
plot(fJMapprox(0:0.1:2),M0f(0:0.1:2))
hold on
scatter(JM,odemax(:,11),50,c,'v','filled')
ylabel('Cell mass (aa)')
xlabel('Intake flux (aa/h)')
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
saveas(gca,fullfile(yourFolder,'Intakeflux_gr_mass'),'png')
%for matlab at home 
%solution based on https://stackoverflow.com/questions/10985946/how-to-export-the-figure-to-color-eps-in-matlab
saveas(gca,fullfile(yourFolder,'Intakeflux_gr_mass.eps'),'epsc')


figure
plot(M0f(0:0.1:2),fJMapprox(0:0.1:2))
hold on
scatter(odemax(:,11),JM,50,c,'filled')
xlabel('Cell mass (aa)')
ylabel('Intake flux (aa/h)')
saveas(gca,fullfile(yourFolder,'Intakeflux_mass'),'png')
saveas(gca,fullfile(yourFolder,'Intakeflux_mass.eps'),'epsc')

figure
% plot(odemax(:,4),1./odemax(:,11),'sq-')
plot(0.4:0.01:2,1./M0f(0.4:0.01:2),'k')
hold on
scatter(odemax(:,4),odemax(:,4)./JM,50,c,'filled')
scatter(gr,gr./(Rpercell.*kselect),50,'kx','LineWidth',1.25)
scatter(dataextract(:,3),dataextract(:,3)./(dataextract(:,1).*dataextract(:,4)*3600),50,'ksq','LineWidth',1.25)
% xlabel('Growth rate (1/h)')
xlabel('\lambda^* (1/h)')
ylabel('\lambda^*/J_M^* (1/aa)')
% ylabel('Growth rate/Intake flux (1/aa)')
% legend('1/M','\lambda/J_M')
xlim([0 2.1])
% ylim([2e-10 20e-10])
ylim([2e-10 16e-10])
% yticks([0.2:0.2:1.6]*1e-9)
xticks(0:0.4:2)
% legend('Theory line','Numerical results','Bremer-Dennis','Hwa')
set(gca,'FontName','Helvetica','FontSize',16)
box on
% export_fig(gca,fullfile(yourFolder,'grJMratio_gr'),'-pdf','-eps')
saveas(gca,fullfile(yourFolder,'grJMratio_gr'),'png')
%for matlab at home 
%solution based on https://stackoverflow.com/questions/10985946/how-to-export-the-figure-to-color-eps-in-matlab
saveas(gca,fullfile(yourFolder,'grJMratio_gr.eps'),'epsc')

%combine with metabolic flux exp data
ft = fittype( 'a*x+b', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [1.96 -40]; %1.96 is exponent from paper www.pnas.org/cgi/doi/10.1073/pnas.1007783107
opts.Robust = 'Bisquare';
opts.StartPoint = [1.96 -10];
[mfit, gof] = fit(log10(mW_aa(:,1)),log10(mW_aa(:,2)),ft, opts )

figure
plot(M0f(0.4:0.1:2),fJMapprox(0.4:0.1:2),'k')
hold on
plot(10.^(8:0.1:10),10.^(mfit(8:0.1:10)),'k--') %[0.49 0.18 0.56]
% scatter(Mpercell*0.7,Rpercell.*kselect,50,'kx') %Bremer - Dennis 
% scatter(dataextract(:,5)*1e-12/maa,dataextract(:,1).*dataextract(:,4)*3600,50,'ksq') 
scatter(mW_aa(:,1),mW_aa(:,2),20,'^','MarkerFaceAlpha',0.75,'MarkerEdgeColor','k')  %[0.64 0.08 0.18]
% hold on 
scatter(odemax(:,11),JM,50,c,'filled')
set(gca,'XScale','log','YScale','log')
xlabel('Log_{10} of cell mass')
ylabel('Log_{10} of intake flux/cell')
legend('DeLong et al.','Numerical results','Location','NorthWest')

% legend('Theory','Experimental fitting','DeLong et al.','Numerical results','Location','NorthWest')
% saveas(gca,fullfile(yourFolder,'flux_mass_scaling'),'png')
% saveas(gca,fullfile(yourFolder,'flux_mass_scaling.eps'),'epsc')
set(gca,'FontName','Helvetica')
box on
ticks = get(gca,'XTickLabel')
expression = '\d*\^\{(\-?\d*)\}'; % Dynamic Regular Expressions
replace = '$1'; 
exponents_X = regexprep(ticks,expression,replace);
ticks = get(gca,'YTickLabel')
exponents_Y = regexprep(ticks,expression,replace);
set(gca,'XTickLabel',exponents_X,'YTickLabel',exponents_Y)
saveas(gca,fullfile(yourFolder,'flux_mass_scaling.pdf'),'pdf')


figure('Position',[10 10 720 510])
plot(M0f(0.3:0.1:2.5),fJMapprox(0.3:0.1:2.5),'k','LineWidth',2)
hold on
plot(10.^(8:0.1:10),10.^(mfit(8:0.1:10)),'k--','LineWidth',2) %[0.49 0.18 0.56]
% scatter(Mpercell*0.7,Rpercell.*kselect,50,'kx') %Bremer - Dennis 
% scatter(dataextract(:,5)*1e-12/maa,dataextract(:,1).*dataextract(:,4)*3600,50,'ksq') 
scatter(mW_aa(:,1),mW_aa(:,2),200,'^','MarkerFaceAlpha',0.75,'MarkerEdgeColor','k')  %[0.64 0.08 0.18]
% hold on 
scatter(odemax(:,11),JM,400,c,'filled')
set(gca,'XScale','log','YScale','log')
xlabel('Log_{10} of cell mass')
ylabel('Log_{10} of intake flux/cell')
legend('DeLong et al.','Numerical results','Location','NorthWest')
box on
set(gca,'FontSize',20)
ticks = get(gca,'XTickLabel')
expression = '\d*\^\{(\-?\d*)\}'; % Dynamic Regular Expressions
replace = '$1'; 
exponents_X = regexprep(ticks,expression,replace);
ticks = get(gca,'YTickLabel')
exponents_Y = regexprep(ticks,expression,replace);
set(gca,'XTickLabel',exponents_X,'YTickLabel',exponents_Y)
set(gca,'FontName','Helvetica')
saveas(gca,fullfile(yourFolder,'flux_mass_scaling_inset.pdf'),'pdf')

figure
plot(M0f(0.2:0.1:2),fJMapprox(0.2:0.1:2),'k')
hold on
scatter(odemax(:,11),JM,50,c,'filled')
scatter(Mpercell*0.6,Rpercell.*kselect,50,'kx','LineWidth',1.25) %Bremer - Dennis 
scatter(dataextract(:,5)*0.6*1e-12/maa,dataextract(:,1).*dataextract(:,4)*3600,50,'k^','LineWidth',1.25) %Hwa-Liu
xlabel('M^* (aa)')
ylabel('J_M^* (aa/h)')
xlim([0 7]*1e9)
xticks([0:1:7]*1e9)
ylim([0 8e9])
yticks([0:8]*1e9)
set(gca,'FontName','Helvetica','FontSize',16)
% ylim([1e8 1e11])
% set(gca,'XScale','log','YScale','log')
% xlabel('Log_{10} of cell mass')
% ylabel('Log_{10} of intake flux/cell')
% legend('Theory line','Numerical results','Bremer - Dennis','Hwa-Liu','Location','NorthWest')
% set(gca,'FontName','Helvetica')
% box on
% ticks = get(gca,'XTickLabel')
% expression = '\d*\^\{(\-?\d*)\}'; % Dynamic Regular Expressions
% replace = '$1'; 
% exponents_X = regexprep(ticks,expression,replace);
% ticks = get(gca,'YTickLabel')
% exponents_Y = regexprep(ticks,expression,replace);
% set(gca,'XTickLabel',exponents_X,'YTickLabel',exponents_Y)
saveas(gca,fullfile(yourFolder,'flux_mass_nutrient_scaling.pdf'),'pdf')
saveas(gca,fullfile(yourFolder,'flux_mass_nutrient_scaling.png'),'png')

figure
scatter(odemax(:,11),JM,'filled')
hold on
plot(M0f(0.2:0.1:2.5),fJMapprox(0.2:0.1:2.5))
xlabel('Cell mass (aa)')
ylabel('Intake flux (aa/h)')
legend('ODEs','Analytical fitting')
set(gca,'XScale','log','YScale','log')
xlim([1e8 1e10])
ylim([1e8 1e11])
ticks = get(gca,'XTickLabel')
expression = '\d*\^\{(\-?\d*)\}'; % Dynamic Regular Expressions
replace = '$1'; 
exponents_X = regexprep(ticks,expression,replace);
ticks = get(gca,'YTickLabel')
exponents_Y = regexprep(ticks,expression,replace);
set(gca,'XTickLabel',exponents_X,'YTickLabel',exponents_Y)
saveas(gca,fullfile(yourFolder,'Intakeflux_mass_loglog'),'png')

gfit = fit(log(M(odemax)),log(JM),'poly1')
% plot(gfit,log(M(maxlist)),log(JM))
figure
scatter(M(odemax),JM,'filled')
hold on
Mlist = logspace(8.8,10,100);
plot(Mlist, exp(gfit.p2)*Mlist.^(gfit.p1))
set(gca,'XScale','log','YScale','log')
xlabel('Cell mass (aa)')
ylabel('Input flux (aa/h)')
saveas(gca,fullfile(yourFolder,'Intakeflux_mass_loglog2'),'png')

figure
bar(odemax(:,4), [odemax(:,5) odemax(:,6)],'stacked')
ylabel('Synthesis allocation')
xlabel('Growth rate (1/h)')
legend('R','P')
saveas(gca,fullfile(yourFolder,'uRuP_gr'),'png')

figure
h = bar(odemax(:,4), odemax(:,7:9),'stacked')
ylabel('Mass fraction')
xlabel('Growth rate (1/h)')
legend('R','P','aa')
h(3).FaceColor = [0.4940 0.1840 0.5560];
h(2).FaceColor = [1 1 0];
saveas(gca,fullfile(yourFolder,'massfr_gr'),'png')

figure
bar(odemax(:,4), [odemax(:,5).*JS odemax(:,6).*JS],'stacked')
hold on
scatter(odemax(:,4), JM,'filled')
ylabel('Synthesis flux (aa/h)')
xlabel('Growth rate (1/h)')
legend('R','P','J_M')
saveas(gca,fullfile(yourFolder,'synthesisflux_gr'),'png')
%for matlab at home 
%solution based on https://stackoverflow.com/questions/10985946/how-to-export-the-figure-to-color-eps-in-matlab
saveas(gca,fullfile(yourFolder,'synthesisflux_gr.eps'),'epsc')

namelist = {'R','P','aa'};
masslist = [NR NP 1];

%theoretical estimation from concentration 
%assume R*ksynmax*[T3]/(KM+[T3]) = S*k2*[aa]/(K+[aa])
% ksynmax_change = S/R*k2; 
% aafr_theoryapprox = @(x) NA/(rho*1e15)*(NR*Ka*x*Ke^2/(2*ksynmax_change*NP)).^(1/3); %x is kmeteff
% gr_theory = @(x) (aafr_theoryapprox(x)*rho*1e15/NA)./(3*NR*Ka/(2*ksynmax_change)+(aafr_theoryapprox(x)*rho*1e15/NA).*(NP./x+NR/ksynmax_change));
% % P_theory = @(x) gr_theory(x).*M0f(gr_theory(x)).*(1+(aafr_theory(x)*rho*1e15/(NA*Ke)).^2)./x;
% aa_theory = @(x) rho*1e15/NA*aafr_theoryapprox(x).*M0f(gr_theory(x));
% P_theory = @(x) (M0f(gr_theory(x))-aa_theory(x)-NR*fR(gr_theory(x)))/NP;
% P_theoryapprox = @(y) (Mx*y+My - NR*fR(y))/NP; %y is growth rate 
% 
% kmetlist = 0:1:0.3*3600;
% aac_exact = NaN*ones(length(kmetlist),1); %[aa]
% gr_exact = NaN*ones(length(kmetlist),1);  %gr*
% fgr = @(x,y) (-NA/(rho*1e15)*x.^2+x)./(NP/(y*Ke^2)*x.^3-NR/ksynmax_change*x.^2+(NR/ksynmax_change+NP/y)*x+NR*Ka/(ksynmax_change));
% syms x
% for i = 1 : length(kmetlist)
%     eqn = 2*NP/(kmetlist(i)*Ke^2)*x^3-NR/ksynmax_change*x^2 - NR*Ka/ksynmax_change == 0;
%     sol = solve(eqn,x,'Real',true);
%     if isempty(sol) ~= 1
%         aac_exact(i,1) = double(sol(double(real(sol))>0));
%         gr_exact(i,1) = fgr(aac_exact(i,1),kmetlist(i));
%     end
% end

%theoretical estimation from concentration 
%assume [T3] = 2*k2*S/R*KM.*[aa]./(ksynmaxT3*Ka)
ksynmax_change = ksynmax; 
Ka_change = KMeff/(2*k2*S/R*KMeff/(ksynmax*Ka));
%x is kmeteff
aafr_theoryapprox = @(x) NA/(rho*1e15)*(NR*Ka_change*x*Ke^2/(2*ksynmax_change*NP)).^(1/3); %x is kmeteff
gr_theory = @(x) (aafr_theoryapprox(x)*rho*1e15/NA)./(3*NR*Ka_change/(2*ksynmax_change)+(aafr_theoryapprox(x)*rho*1e15/NA).*(NP./x+NR/ksynmax_change));
aa_theory = @(x) rho*1e15/NA*aafr_theoryapprox(x).*M0f(gr_theory(x));
% P_theory = @(x) gr_theory(x).*M0f(gr_theory(x)).*(1+(aafr_theory(x)*rho*1e15/(NA*Ke)).^2)./x;
P_theory = @(x) (M0f(gr_theory(x))-aa_theory(x)-NR*fR(gr_theory(x)))/NP;
% P_theoryapprox = @(y) (Mx*y+My - NR*fR(y))/NP; %y is growth rate 

for i = 1 : length(namelist)
    yL = sprintf('%s^* (molecules/cell)',namelist{i}); 
    yR = ['\phi_' sprintf('{%s}^*',namelist{i})]; 
    figure
    yyaxis left
%     plot(odemax(:,4),odemax(:,i),'o-','MarkerFaceColor',[0 0.4470 0.7410])
%     plot(odemax(:,4),odemax(:,i),'o-','MarkerFaceColor','w')
    switch i
        case 1
            plot(0:0.1:2, fR(0:0.1:2),'LineWidth',2,'Color','k')
            yticks(0:2e4:10e4)
            ylim([0 10e4])
        case 2
            plot(gr_theory(0:0.1:0.3*3600),P_theory(0:0.1:0.3*3600),'Color','k')
            ylim([(floor(min(odemax(:,2))/1e6))*0.1*1e6 (round(max(odemax(:,2))/1e6)+1)*1e6])
        case 3
            plot(gr_theory(0:0.1:0.3*3600),aa_theory(0:0.1:0.3*3600),'Color','k')
            ylim([(floor(min(odemax(:,3))/1e5))*0.5*1e5 (round(max(odemax(:,3))/1e6)+0.5)*1e6])
    end
    hold on
    if i == 1
        scatter(odemax(:,4),odemax(:,i),100,c,'filled')
    else
        scatter(odemax(:,4),odemax(:,i),50,c,'filled')
    end
    ylabel(yL)
    set(gca,'FontName','Helvetica')
    yyaxis right
    switch i
        case 1
            plot(0:0.1:2, hfit(0:0.1:2),'LineWidth',2)
            ylim([0 0.4])
        case 2
            plot(gr_theory(0:0.1:0.3*3600),NP*P_theory(0:0.1:0.3*3600)./M0f(gr_theory(0:0.1:0.3*3600)))
            ylim([0.2 1])
        case 3
%             plot(gr_exact,aac_exact*NA/(rho*1e15))    
%             plot(gr_theory(0:0.1:0.3*3600),aafr_theoryapprox(0:0.1:0.3*3600))
            plot(gr_theory(0:0.1:0.3*3600),1-hfit(gr_theory(0:0.1:0.3*3600))-NP*P_theory(0:0.1:0.3*3600)'./M0f(gr_theory(0:0.1:0.3*3600))')
    end
    hold on
    if i == 1
        scatter(odemax(:,4),mfraction(:,i),100,c,'v','filled')
    else
        scatter(odemax(:,4),mfraction(:,i),50,c,'v','filled')
    end
    xlabel('\lambda^* (1/h)')
    xlim([0 2])
    xticks(0:0.2:2)

    ylabel(yR) 
    if i == 3
        ylim([floor(min(odemax(:,9)*1e4))/2*1e-4 round(max(odemax(:,9)*1e3))*1e-3])
    end
    if i == 1        
        xticks(0:0.4:2)
        set(gca,'FontSize',20)
    end
    set(gca,'FontName','Helvetica')
    ax = gca;
    ax.YAxis(1).Color = 'k';
%     ax.YAxis(2).Color = 'k';
    saveas(gca,fullfile(yourFolder,namelist{i}),'png')
    %for matlab at home 
    %solution based on https://stackoverflow.com/questions/10985946/how-to-export-the-figure-to-color-eps-in-matlab
    saveas(gca,fullfile(yourFolder,[namelist{i} '.eps']),'epsc')

end

%compare with analytical calculated from concentration approximation
%if needed then can restore and check 

% kmetlist = 0:1:0.3*3600;
% aac_exact = NaN*ones(length(kmetlist),1); %[aa]
% gr_exact = NaN*ones(length(kmetlist),1);  %gr*
% R_exact = NaN*ones(length(kmetlist),1);  %[R]
% P_exact = NaN*ones(length(kmetlist),1);  %[P]
% % x - [aa], y - kmet 
% fgr = @(x,y) (-NA/(rho*1e15)*x.^2+x)./(NP/(y*Ke^2)*x.^3-NR/ksynmax_change*x.^2+(NR/ksynmax_change+NP/y)*x+NR*Ka_change/(ksynmax_change));
% Pc_f = @(x,y) (fgr(x,y).*rho*1e15./NA)./(y./(1+(x/Ke).^2)); %[P] theory 
% Rc_f = @(x,y) (rho*1e15/NA-x-Pc_f(x,y)*NP)/NR;
% 
% syms x
% for i = 1 : length(kmetlist)    
%     eqn = 2*NP/(kmetlist(i)*Ke^2)*x^3-(NR/ksynmax_change)/(rho*1e15/NA)*x^2 - NR*Ka_change/ksynmax_change == 0;
%     sol = solve(eqn,x,'Real',true);
%     if isempty(sol) ~= 1
%         aac_exact(i,1) = double(sol(double(real(sol))>0));
%         gr_exact(i,1) = fgr(aac_exact(i,1),kmetlist(i));
%         P_exact(i,1) = Pc_f(aac_exact(i,1),kmetlist(i));
%         R_exact(i,1) = Rc_f(aac_exact(i,1),kmetlist(i));
%     end
% end
% 
% for i = 1 : length(namelist)
%     yL = sprintf('#%s/cell',namelist{i}); 
%     yR = sprintf('%s mass fraction',namelist{i}); 
%     figure
%     yyaxis left
%     switch i
%         case 1
%             plot(gr_exact,M0f(gr_exact).*R_exact*NA/(rho*1e15),'LineWidth',2,'Color','k')
%             yticks(0:2e4:10e4)
%             ylim([0 10e4])
%         case 2
%             plot(gr_exact, M0f(gr_exact).*P_exact*NA/(rho*1e15),'LineWidth',2)
%             ylim([(floor(min(odemax(:,2))/1e6))*0.1*1e6 (round(max(odemax(:,2))/1e6)+1)*1e6])
%         case 3
%             plot(gr_exact,M0f(gr_exact).*aac_exact*NA/(rho*1e15))    
%             ylim([(floor(min(odemax(:,3))/1e5))*0.5*1e5 (round(max(odemax(:,3))/1e6)+1)*1e6])
%     end
%     hold on
%     if i == 1
%         scatter(odemax(:,4),odemax(:,i),100,c,'filled')
%     else
%         scatter(odemax(:,4),odemax(:,i),50,c,'filled')
%     end
%     ylabel(yL)
%     set(gca,'FontName','Helvetica')
%     yyaxis right
% %     plot(odemax(:,4),mfraction(:,i),'sq-','MarkerFaceColor',[0.8500 0.3250 0.0980])
% %     plot(odemax(:,4),mfraction(:,i),'--','MarkerFaceColor','w')
%     switch i
%         case 1
%             plot(gr_exact, R_exact*NR*NA/(rho*1e15),'LineWidth',2)
%             ylim([0 0.4])
%         case 2
%             plot(gr_exact, P_exact*NP*NA/(rho*1e15),'LineWidth',2)
%             ylim([0.1 1])
%         case 3
%             plot(gr_exact,aac_exact*NA/(rho*1e15))    
%     end
%     hold on
%     if i == 1
%         scatter(odemax(:,4),mfraction(:,i),100,c,'v','filled')
%     else
%         scatter(odemax(:,4),mfraction(:,i),50,c,'v','filled')
%     end
%     xlabel('Growth rate (1/h)')
%     xlim([0 2])
%     xticks(0:0.2:2)
% 
%     ylabel(yR) 
%     if i == 3
%         ylim([floor(min(odemax(:,9)*1e4))/2*1e-4 round(max(odemax(:,9)*1e3))*2e-3])
%     end
%     if i == 1        
%         xticks(0:0.4:2)
%         set(gca,'FontSize',20)
%     end
%     set(gca,'FontName','Helvetica')
%     ax = gca;
%     ax.YAxis(1).Color = 'k';
% %     ax.YAxis(2).Color = 'k';
%     saveas(gca,fullfile(yourFolder,[namelist{i} 'analytical_vs_ode']),'png')
%     %for matlab at home 
%     %solution based on https://stackoverflow.com/questions/10985946/how-to-export-the-figure-to-color-eps-in-matlab
%     saveas(gca,fullfile(yourFolder,[namelist{i} 'analytical_vs_ode.eps']),'epsc')
% 
% end


figure
scatter(odemax(:,4),kmetsample*1/3600,'filled')
xlabel('Growth rate (1/h)')
ylabel('k^{eff}_{met} (1/s)')
saveas(gca,fullfile(yourFolder,'kmeteff_gr'),'png')

%for all values of kmetsample and grsample
% figure
% scatter(kmetsample.*0.05./(0.4*3600-kmetsample),odemax(:,4),'filled')
% ylabel('Growth rate (1/h)')
% xlabel('[Nutrients] (mM)')

figure
scatter(kmetsample.*0.4./(0.3*3600-kmetsample),odemax(:,4),'d','filled')
ylabel('Growth rate (1/h)')
xlabel('[Nutrients] (mM)')
saveas(gca,fullfile(yourFolder,'kmeteff_nutrient'),'png')

%fit kmet* vs gr*
% x0 = [0.2*3600 0.01];
% fitfun = fittype(@(a,b,x) a*x./(b+x) );
% [fitted_curve,gof] = fit(kmetsample.*0.05./(0.4*3600-kmetsample),kmetsample,fitfun,'StartPoint',x0);
% coeffvals_kmet = coeffvalues(fitted_curve);
% figure
% scatter(kmetsample.*0.05./(0.4*3600-kmetsample),kmetsample/3600,'filled')
% hold on
% plot(0:0.01:0.2,fitted_curve(0:0.01:0.2)/3600)

x0 = [0.2*3600 0.01];
fitfun = fittype(@(a,b,x) a*x./(b+x) );
[fitted_curve,gof] = fit(kmetsample.*0.4./(0.3*3600-kmetsample),kmetsample,fitfun,'StartPoint',x0);
coeffvals_kmet = coeffvalues(fitted_curve);

% c = [0 0 0; 96 60 38; 142 88 56; 193 151 96; 255 199 127]/255; %list of color
c = [77 191 237; 113 152 241; 148 114 244; 188 71 248; 255 0 255]/255; %list of color

% figure('Position',[1 1 587 521])
figure
plot(0:0.01:2,fitted_curve(0:0.01:2)/3600,'k')
hold on
scatter(kmetsample.*0.4./(0.3*3600-kmetsample),kmetsample/3600,50,c,'filled')
ylabel('k^{eff}_{met} (1/s)')
xlabel('Nutrient quality')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
xticks([0 2])
xticklabels({'Poor','Rich'})
% colormap('copper')
% colorbar('southoutside','TickLabels',[],'Ticks',[])
saveas(gca,fullfile(yourFolder,'kmeteff_quality'),'png')
%for matlab at home 
%solution based on https://stackoverflow.com/questions/10985946/how-to-export-the-figure-to-color-eps-in-matlab
saveas(gca,fullfile(yourFolder,'kmeteff_quality.eps'),'epsc')

close all
end