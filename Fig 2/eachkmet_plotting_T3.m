 function eachkmet_plotting_T3(odestore,para,j,V,M,kmetsample,Foldername,ktrans_T3)
Ka = para(1);
Ke = para(2);
NA = para(3);
ksynmax = para(4);
kmet = kmetsample(j);

%illustrate flux limited by transporter or ribosomes
% fpath = 'D:\Onedrive-NTU\OneDrive - Nanyang Technological University\PhD\Kinetics\Kinetics methods\OPTI-master\Figures\growth_HB\';
% yourFolder = [fpath,sprintf('kmet = %.2f',kmet/3600)];
% yourFolder = ['../OPTI-master/Figures/growth_HB/' sprintf('kmet = %.2f',kmet/3600)];
yourFolder = ['../growth_HB/' Foldername];

if exist(yourFolder, 'dir') ~= 7 %folder does not exist
       mkdir(yourFolder)
end

figure
JMstore = odestore{j}(:,2).*kmetsample(j)./(1+(odestore{j}(:,3)./(Ke*V(odestore{j})*1e-15*NA)).^2);
JSstore =  odestore{j}(:,1).*ktrans_T3(odestore{j});
JMstore./(JSstore+odestore{j}(:,4).*odestore{j}(:,3)); %check flux balance
[val, id] = max(odestore{j}(:,4)); %max growth rate
frange = find(abs(JMstore-JMstore(id))<=0.3*JMstore(id)); %find within certain range
plot(odestore{j}(frange,7),[JMstore(frange) JSstore(frange)])
hold on
plot(odestore{j}(frange,7), [odestore{j}(frange,2).*kmetsample(j) odestore{j}(frange,1).*ksynmax],'--')
title(sprintf('k{met}^{eff} = %.2f 1/s',kmetsample(j)/3600))
xlabel('R mass fraction')
ylabel('Flux (aa/h)')
legend('J_M','J_S','J_M^{max}','J_S^{max}','Location','northwest')
saveas(gca,fullfile(yourFolder,'flux drives'),'png')

% xlim([0.1 0.25])
% ylim([1e9 2e9])

figure
plot(odestore{j}(frange,7), odestore{j}(frange,4)./JMstore(frange))
hold on
plot(odestore{j}(frange,7), 1./M(odestore{j}(frange,:)))
xlabel('R mass fraction')
ylabel('\lambda/J_M (1/aa)')
legend('\lambda/J_M','1/M')
ylim([1e-10 1e-9])
saveas(gca,fullfile(yourFolder,'gr_JM ratio'),'png')

%pick to illustrate which partition to choose via flux
[val, id] = max(odestore{j}(:,4)); %max growth rate
flux_pick = NaN*ones(3,6);
%1st - 2 - 3 - 4 - 5 col: R mass fraction - uR - JM - JS - JMmax - JSmax
rgap = [-5 0 5];
for i = 1 : 3
    if i == 3 
        if id - rgap(i) >= 1
            flux_pick(i,:) =  [odestore{j}(id-rgap(i),7) odestore{j}(id-rgap(i),5) JMstore(id-rgap(i)) JSstore(id-rgap(i)) odestore{j}(id-rgap(i),2).*kmetsample(j) odestore{j}(id-rgap(i),1).*ksynmax]; 
        else       
            flux_pick(i,:) =  [odestore{j}(id-rgap(i)+1,7) odestore{j}(id-rgap(i)+1,5) JMstore(id-rgap(i)+1) JSstore(id-rgap(i)+1) odestore{j}(id-rgap(i)+1,2).*kmetsample(j) odestore{j}(id-rgap(i)+1,1).*ksynmax]; 
        end
    else
        flux_pick(i,:) =  [odestore{j}(id-rgap(i),7) odestore{j}(id-rgap(i),5) JMstore(id-rgap(i)) JSstore(id-rgap(i)) odestore{j}(id-rgap(i),2).*kmetsample(j) odestore{j}(id-rgap(i),1).*ksynmax]; 
    end
end

figure
b = bar(flux_pick(:,1),[flux_pick(:,2).*flux_pick(:,4) (1-flux_pick(:,2)).*flux_pick(:,4)],'stacked')
b(1).FaceColor = [0 0 139]/225;
b(2).FaceColor = [1 1 1];
hold on
scatter(flux_pick(:,1),flux_pick(:,3),'filled')
plot(flux_pick(:,1),flux_pick(:,5:6),'sq--')
legend('R','P','J_M','J_M^{max}','J_S^{max}','Location','northwest')
xlabel('R mass fraction')
ylabel('Flux (aa/h)')
saveas(gca,fullfile(yourFolder,'flux_mass fraction'),'png')
saveas(gca,fullfile(yourFolder,'flux_mass fraction.pdf'),'pdf')


%and same mass with different ribosome fraction
figure
store = odestore{j};
h = area(store(frange,7), [store(frange,7).*M(store(frange,:)) store(frange,8).*M(store(frange,:)) store(frange,9).*M(store(frange,:))])
xlabel('R mass fraction')
ylabel('Cell mass (aa)')
legend('R','P','aa')
h(1).FaceColor = [0 0 139]/225;
h(2).FaceColor = [1 1 1];
h(3).FaceColor = [0.4940 0.1840 0.5560];
saveas(gca,fullfile(yourFolder,'cellmass_partition'),'png')
saveas(gca,fullfile(yourFolder,'cellmass_partition.pdf'),'pdf')

figure
store = odestore{j};
h = area(store(frange,7), [store(frange,7).*M(store(frange,:)) store(frange,8).*M(store(frange,:)) store(frange,9).*M(store(frange,:))])
xlabel('R mass fraction')
ylabel('Cell mass (aa)')
h(1).FaceColor = [0 0 139]/225;
h(2).FaceColor = [1 1 1];
h(3).FaceColor = [0.4940 0.1840 0.5560];
ylim([2.26 2.27]*1e9)
xlim([0.1 0.55])
set(gca,'FontSize',25)

% xlim([0.1 0.25])
%and growth rate with different ribosome fraction
figure
yyaxis left
plot(odestore{j}(frange,7),odestore{j}(frange,4))
ylabel('Growth rate (1/h)')
yyaxis right
plot(odestore{j}(frange,7),log(2)./odestore{j}(frange,4))
xlabel('R mass fraction')
ylabel('Doubling time (h)')
saveas(gca,fullfile(yourFolder,'phiR_gr'),'png')

% xlim([0.1 0.25])
%synthesis allocation vs ribosome mass fraction
figure
area(odestore{j}(frange,7),[odestore{j}(frange,5) odestore{j}(frange,6)])
xlabel('R mass fraction')
ylabel('Synthesis allocation')
legend('R','P')
saveas(gca,fullfile(yourFolder,'phiR_partition'),'png')

% xlim([0.1 0.25])
%#R,P,aa vs ribosome mass fraction
% semilogy(odestore{j}(frange,7),odestore{j}(frange,1:3))
namelist = {'R','P','aa'};
for i = 1 : 3
figure
plot(odestore{j}(frange,7),odestore{j}(frange,i))
xlabel('R mass fraction')
ylabel(sprintf('#%s (molecules)',namelist{i}))
saveas(gca,fullfile(yourFolder,namelist{i}),'png')

% xlim([0.1 0.25])
end
close all