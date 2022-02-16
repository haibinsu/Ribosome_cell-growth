function fluxdrive_allkmetT3(odestore,para,V,M,kmetsample,odemax,nstore,Foldername,ktrans_T3,glawdata)
Ka = para(1);
Ke = para(2);
NA = para(3);
ksynmax = para(4);

%illustrate flux limited by transporter or ribosomes
% fpath = 'D:\Onedrive-NTU\OneDrive - Nanyang Technological University\PhD\Kinetics\Kinetics methods\OPTI-master\Figures\growth_HB\';
% yourFolder = [fpath,sprintf('kmet = %.2f',kmet/3600)];
% yourFolder = ['../OPTI-master/Figures/growth_HB/' 'allkmetsameplot'];
yourFolder = ['../growth_HB/', Foldername];
if exist(yourFolder, 'dir') ~= 7 %folder does not exist
       mkdir(yourFolder)
end

figure
% c = [0 0 0; 96 60 38; 142 88 56; 193 151 96; 255 199 127]/255; %list of color
c = [77 191 237; 113 152 241; 148 114 244; 188 71 248; 255 0 255]/255; %list of color
for j =  1 : length(kmetsample)
    JMstore = odestore{j}(:,2).*kmetsample(j)./(1+(odestore{j}(:,3)./(Ke*V(odestore{j})*1e-15*NA)).^2);
    [val, id] = max(odestore{j}(:,4)); %max growth rate
    frange = find(abs(JMstore-JMstore(id))<=0.1*JMstore(id));
%     frange = find(JMstore-JMstore(id) <= -0.01*JMstore(id) | JMstore-JMstore(id) <= 0.01*JMstore(id));
%     if j == 1 | j == 2 | j == 3
%         frange = find( val-odestore{j}(:,4) >= 0.1*val & odestore{j}(:,5) <= 0.2);
%         frange = find(JMstore-JMstore(id) <= -0.1*JMstore(id) & odestore{j}(:,5) <= 0.2);

%     end
    %     plot(JMstore, odestore{j}(:,4))
%     plot(odestore{j}(frange,7), JMstore(frange))
    plot(odestore{j}(frange,7),odestore{j}(frange,4),'Color',c(j,:))
    hold on
end
hfit = fit(odemax(:,7),odemax(:,4),'poly1')
plot((0.1:0.1:0.5),hfit.p1*(0.1:0.1:0.5)+hfit.p2,'Color','k')
scatter(odemax(:,7),odemax(:,4),50,c,'filled')
% hold on
% scatter(glawdata(:,2)./min(glawdata(:,2))*min(odemax(:,7)),glawdata(:,1))
% glawfit = fit(glawdata(:,2),glawdata(:,1),'poly1')
% plot((0.12:0.1:0.5),glawfit.p1*(0.12:0.1:0.5)+glawfit.p2,'k--')
scatter(glawdata(1:6,2),glawdata(1:6,1),'k','diamond','LineWidth',1.25) %Hwa
scatter(glawdata(7:11,2),glawdata(7:11,1),50,'k','<','LineWidth',1.25)  %Forchhammer & Lindahl
scatter(glawdata(12:end,2),glawdata(12:end,1),'k','o','LineWidth',1.25) %Bremer & Dennis
ylabel('\lambda^* (1/h)')
xlabel('\phi_R^*')
xlim([0.1 0.5])
ylim([0 2])
yticks([0:0.4:2])
xticks([0.1:0.1:0.5])
set(gca,'FontName','Arial','FontSize',16)
saveas(gca,fullfile(yourFolder,'gr_phiR'),'png')
saveas(gca,fullfile(yourFolder,'phiR_gr'),'fig')
% saveas(gca,fullfile(yourFolder,'phiR_gr'),'eps')
%for matlab at home 
%solution based on https://stackoverflow.com/questions/10985946/how-to-export-the-figure-to-color-eps-in-matlab
saveas(gca,fullfile(yourFolder,'phiR_gr.eps'),'epsc')


figure
scatter(odemax(:,2).*kmetsample./(1+(odemax(:,3)./(Ke*V(odemax)*1e-15*NA)).^2),odemax(:,4),[],c,'filled')
xlabel('Intake flux (aa/h)')
ylabel('Growth rate (1/h)')
box on
saveas(gca,fullfile(yourFolder,'flux_gr_only'),'png')
saveas(gca,fullfile(yourFolder,'flux_gr_only'),'fig')
%for matlab at home 
%solution based on https://stackoverflow.com/questions/10985946/how-to-export-the-figure-to-color-eps-in-matlab
saveas(gca,fullfile(yourFolder,'flux_gr_only.eps'),'epsc')


%to illustrate input fluxes drive growth
% c = [0 0 0; 96 60 38; 142 88 56; 193 151 96; 255 199 127]/255; %list of color
c = [77 191 237; 113 152 241; 148 114 244; 188 71 248; 255 0 255]/255; %list of color

figure
for j =  1 : length(kmetsample)
    JMstore = odestore{j}(:,2).*kmetsample(j)./(1+(odestore{j}(:,3)./(Ke*V(odestore{j})*1e-15*NA)).^2);
    [val, id] = max(odestore{j}(:,4)); %max growth rate
    frange = find(abs(JMstore-JMstore(id))<=0.3*JMstore(id));
    plot(JMstore(frange),odestore{j}(frange,4),'Color',c(j,:))
    hold on
end
% c = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.4940 0.1840 0.5560; 0.4660 0.6740 0.1880]; %list of color
% c = 1:5;
s = scatter(odemax(:,2).*kmetsample./(1+(odemax(:,3)./(Ke*V(odemax)*1e-15*NA)).^2),odemax(:,4),[],c,'filled')
% colormap(gca,'copper')
% s.CData
xlabel('Intake flux (aa/h)')
ylabel('Growth rate (1/h)')
saveas(gca,fullfile(yourFolder,'flux_gr'),'png')
saveas(gca,fullfile(yourFolder,'flux_gr'),'fig')

%for matlab at home 
%solution based on https://stackoverflow.com/questions/10985946/how-to-export-the-figure-to-color-eps-in-matlab
saveas(gca,fullfile(yourFolder,'flux_gr.eps'),'epsc')

figure
% colorscheme = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880]};
yyaxis left
for j =  1 : length(kmetsample)
    JMstore = odestore{j}(:,2).*kmetsample(j)./(1+(odestore{j}(:,3)./(Ke*V(odestore{j})*1e-15*NA)).^2);
    [val, id] = max(odestore{j}(:,4)); %max growth rate
    frange = find(abs(JMstore-JMstore(id))<=0.05*JMstore(id));
    plot(JMstore(frange),odestore{j}(frange,4),'-','Color',c(j,:))
    hold on
end
% c = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.4940 0.1840 0.5560; 0.4660 0.6740 0.1880]; %list of color
scatter(odemax(:,2).*kmetsample./(1+(odemax(:,3)./(Ke*V(odemax)*1e-15*NA)).^2),odemax(:,4),[],c,'filled')
ylabel('Growth rate (1/h)')
yyaxis right
scatter(odemax(:,2).*kmetsample./(1+(odemax(:,3)./(Ke*V(odemax)*1e-15*NA)).^2),kmetsample/3600,[],c,'s','filled')
xlabel('Intake flux (aa/h)')
ylabel('k_{met}^{eff} (1/s)')
ylim([0 0.35])
saveas(gca,fullfile(yourFolder,'flux_gr_kmet'),'png')

figure
scatter(1./odemax(:,11),odemax(:,4)./(odemax(:,2).*kmetsample./(1+(odemax(:,3)./(Ke*V(odemax)*1e-15*NA)).^2)),'filled')
hold on
plot(1./odemax(:,11),1./odemax(:,11),'sq')
plot(10.^(-10:0.1:-9),10.^(-10:0.1:-9))
xlabel('Inverse of cell mass (1/aa)')
ylabel('\lambda/J_M (1/aa)')

figure
plot(odemax(:,4),1./odemax(:,11),'sq-')
hold on
scatter(odemax(:,4),odemax(:,4)./(odemax(:,2).*kmetsample./(1+(odemax(:,3)./(Ke*V(odemax)*1e-15*NA)).^2)),[],c,'filled')
xlabel('Growth rate (1/h)')
ylabel('\lambda/J_M (1/aa)')
legend('1/M','\lambda/J_M')
xlim([0 2])
xticks(0:0.2:2)
box on
saveas(gca,fullfile(yourFolder,'grJMratio_gr'),'png')
%for matlab at home 
%solution based on https://stackoverflow.com/questions/10985946/how-to-export-the-figure-to-color-eps-in-matlab
saveas(gca,fullfile(yourFolder,'grJMratio_gr.eps'),'epsc')

figure
scatter(kmetsample/3600,odemax(:,4)./(odemax(:,2).*kmetsample./(1+(odemax(:,3)./(Ke*V(odemax)*1e-15*NA)).^2)),'filled')
hold on
plot(kmetsample/3600,1./odemax(:,11),'sq-')
xlabel('k_{met}^{eff} (1/s)')
ylabel('\lambda/J_M (1/aa)')
legend('\lambda/J_M','1/M')
saveas(gca,fullfile(yourFolder,'grJMratio_kmet'),'png')

figure
yyaxis left
scatter(odemax(:,7),log(2)./odemax(:,4),'filled')
ylabel('Doubling time (h)')
yyaxis right
scatter(odemax(:,7),odemax(:,4),'sq','filled')
ylabel('Growth rate (1/h)')
xlabel('R mass fraction')
saveas(gca,fullfile(yourFolder,'phiR_dbltime_gr'),'png')

figure
for j =  1 : length(kmetsample)
    JMstore = odestore{j}(:,2).*kmetsample(j)./(1+(odestore{j}(:,3)./(Ke*V(odestore{j})*1e-15*NA)).^2);
    JSstore =  odestore{j}(:,1).*ktrans_T3(odestore{j});

    [val, id] = max(odestore{j}(:,4)); %max growth rate
    frange = find(abs(JMstore-JMstore(id))<=0.3*JMstore(id));
    %     plot(JMstore, odestore{j}(:,4))
%     plot(odestore{j}(frange,7), JMstore(frange))
    plot(odestore{j}(frange,7),log(2)./odestore{j}(frange,4))
    hold on
end
scatter(odemax(:,7),log(2)./odemax(:,4),'filled')
ylabel('Doubling time (h)')
xlabel('R mass fraction')
saveas(gca,fullfile(yourFolder,'phiR_dbltime'),'png')


figure
for j =  1 : length(kmetsample)
    JMstore = odestore{j}(:,2).*kmetsample(j)./(1+(odestore{j}(:,3)./(Ke*V(odestore{j})*1e-15*NA)).^2);
    JSstore =  odestore{j}(:,1).*ktrans_T3(odestore{j});

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
saveas(gca,fullfile(yourFolder,'flux_phiR'),'png')



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
scatter(odemax(:,2).*kmetsample./(1+(odemax(:,3)./(Ke*V(odemax)*1e-15*NA)).^2),odemax(:,11),'filled')
xlabel('Input flux (aa/h)')
ylabel('Cell mass (aa)')
saveas(gca,fullfile(yourFolder,'flux_mass'),'png')

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
   xlabel('R mass fraction')
   saveas(gca,fullfile(yourFolder,sprintf('flux_gr_phiR_kmetid = %d',j)),'png')

end
close all
end