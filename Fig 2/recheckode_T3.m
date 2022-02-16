function recheckode_T3(j,odemax,kmetsample,flist,nstore,para2,M,Rmfc,Pmfc,aamfc,V,ktrans_T3,Foldername,Vf,T3convertf_T3)
ksynmax = para2(1);
nSim = para2(2);
NR = para2(3);
NP = para2(4);
Ke = para2(5);
NA = para2(6);
% yourFolder = ['../OPTI-master/Figures/growth_HB/' sprintf('kmet = %.2f_timeode',kmet/3600)];
yourFolder = ['../growth_HB/' Foldername];

if exist(yourFolder, 'dir') ~= 7 %folder does not exist
       mkdir(yourFolder)
end

%to recheck ode smiluation after par-for
kmet = kmetsample(j);
%metabolic rate
kmetf = @(y) kmet./(1+(y(3)/(Ke*NA*Vf(y)*1e-15))^2);
flist{5} = kmetf; 
store = NaN*ones(1,10);
%1st - 10th col: R - P - aa - gr - uR - uP - phiR - phiP - phiaa - T3

Mth = nstore(j,6); %2*Mth = mass threshold for cell division
R0 = odemax(j,1);
P0 = odemax(j,2);
M0 = odemax(j,11); %aa
Mth/M0
aa0 = M0 - (NR*R0 + P0*NP);
y0 = [R0 P0 aa0];
par = [ksynmax kmet odemax(j,5) NR NP];

k = 1;
t_tot = 0;
y_tot = y0;
tspan = [0 200];
gr = 1;
Opt = odeset('Events', @(t,y) myEvent_growth_opt2(t,y,Mth,par));

ncell = NaN*ones(nSim,1);
ncell(1) = 1;
Mtotal = NaN*ones(nSim,1);
Mtotal(1) = ncell(1)*M0;
time = NaN*ones(nSim,1);
time(1) = 0;    
M0store = NaN*ones(nSim,1);
M0store(1) = M0;

while k <= nSim
    [t,y, te, ye, ie] = ode15s(@(t,y) ode_growth_HB_T3(t,y,par,flist),tspan,y0, Opt);
    if isempty(ye) ~= 1
        t_tot = [t_tot ; [t;te] + t_tot(end)];
        y_tot = [ y_tot ; [y;ye] ];    
        y0 = ye/2;
        k = k + 1;
        ncell(k) = ncell(k-1)*2;
        Mtotal(k) = ncell(k-1)*2*sum(y0.*[NR NP 1]);
        M0store(k) = sum(y0.*[NR NP 1]);
        time(k) = t_tot(end);
    else
        gr = NaN;
        break;
    end
end

%to calculate cell number and total mass in exponential steady state
if isnan(gr) ~= 1 
    [t,y,te,ye,ie] = ode15s(@(t,y) ode_growth_HB_T3(t,y,par,flist),0:1e-3:2*te,y0,Opt);
    semilogy(t,y)
    store = [y0 log(2)/t(end) par(3) 1 - par(3) Rmfc(y0) Pmfc(y0) aamfc(y0) T3convertf_T3(y0)];
end

figure
semilogy(t_tot,y_tot)
xlabel('Time (h)')
ylabel('# molecules')
legend('R','P','aa')
saveas(gca,fullfile(yourFolder,'molecule_time'),'png')

%consistent check 
odemax(j,4)/(log(2)/t(end));

figure
yyaxis left
plot(time,M0store)
ylim([1e8 5e9])
ylabel('Initial cell mass M_0 (aa)')
yyaxis right
semilogy(time,ncell)
ylabel('# cell')
xlabel('Time (h)')

figure
semilogy(time,ncell)
hold on
semilogy(time,Mtotal)
legend('# cell','Total cell mass (aa)')
xlabel('Time (h)')
ylabel('Counts')
saveas(gca,fullfile(yourFolder,'number_mass_time'),'png')


fit(time(20:end), Mtotal(20:end),'exp1')
fit(time(20:end), ncell(20:end),'exp1')
    
%evaluate dy/dt 
dy_dt = NaN*ones(length(t),3);
for i = 1:length(t)
    dy_dt(i,:) = ode_growth_HB_T3(t(i),y(i,:),par,flist);
end

%time dependent input flux 
JM = y(:,2)*kmet./(1+(y(:,3)./(Ke*V(y)*NA*1e-15)).^2);
%time dependent synthesis flux 
JS = y(:,1).*ktrans_T3(y);
% JS./(y(:,1).*ksynmax.*y(:,3)./(Ka*V(y)*NA*1e-15+y(:,3)))

hfit = fit(M(y),JM,'poly1');
hfit.p1/((log(2)/t(end)))

figure
plot(hfit,M(y),JM)
xlabel('M(t) (aa)')
ylabel('J_M(t) (#aa/h)')

JM./(JS+dy_dt(:,3));
JS./(NR*dy_dt(:,1)+NP*dy_dt(:,2));

figure
plot(t,JM)
hold on
plot(t,JS)
plot(t,NR*dy_dt(:,1)+NP*dy_dt(:,2)+dy_dt(:,3))
legend('J_M','J_S','dM/dt')
saveas(gca,fullfile(yourFolder,'time_JM'),'png')

figure
plot(t,[JM,JS])
ylabel('Flux (#aa/h)')
xlabel('Time (h)')
legend('J_M','J_S')

figure
area(t,[par(3)*JS (1-par(3))*JS])
xlabel('Time (h)')
legend('R','P')
ylabel('Synthesis flux J_S (#aa/h)')
saveas(gca,fullfile(yourFolder,'time_JS partition'),'png')

figure
h = area(t,[NR*y(:,1)./M(y), NP*y(:,2)./M(y), y(:,3)./M(y)])
xlabel('Time (h)')
legend('R','P','aa')
ylabel('Mass fraction')
h(3).FaceColor = [0.4940 0.1840 0.5560];


namelist = {'R','P','aa'};
masslist = [NR NP 1];

for i = 1 : length(namelist)
    yL = sprintf('#%s',namelist{i}); 
    yR = sprintf('d%s/dt',namelist{i}); 
    figure
    yyaxis left
    plot(t,y(:,i))
    ylabel(yL)
    yyaxis right
    plot(t,dy_dt(:,i))
    xlabel('Time (h)')
    ylabel(yR) 
    
%     figure
%     plot(t,y(:,i))
%     hold on
%     plot(t,dy_dt(:,i))
%     xlabel('Time (h)')
%     ylabel('# or fluxes (aa/h)') 
%     legend(yL,yR)
end

%check growth rate
for i = 1 : length(masslist)
    plot(t,dy_dt(:,i)./y(:,i))
    hold on
end
xlabel('Time (h)')
ylabel('1/XdX/dt') 
legend('R','P','aa')
ylim([0 2])

namelist = {'R','P'};
ulist = [par(3) 1-par(3)];
allocation = repmat(ulist,length(t),1);

for i = 1 : length(ulist)
    yL = sprintf('%s mass fraction',namelist{i}); 
    yR = sprintf('%s synthesis allocation',namelist{i}); 
    figure
    yyaxis left
    plot(t,masslist(i)*y(:,i)./M(y))
    ylabel(yL)    
    ylim([0 1])

    yyaxis right
    plot(t,allocation(:,i))
    xlabel('Time (h)')
    ylabel(yR) 
    ylim([0 1])
end

figure
area(t,allocation)
xlabel('Time (h)')
legend('R','P')
ylabel('Synthesis allocation')


figure
yyaxis left
plot(t,M(y))
ylabel('Cell mass (#aa)')
yyaxis right
plot(t,V(y))
ylabel('Cell volume (\mum^3)')
xlabel('Time (h)')

close all
end