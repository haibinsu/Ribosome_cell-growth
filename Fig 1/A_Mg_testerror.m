Mgpep = [1.3; 2.3; 3.4; 4.6; 7.5];
kcatKMc = [60; 86; 117; 147; 167; 180];  % cognate AAA uM^-1s^-1
kcatKMnc = [19; 36; 66; 139; 327; 1750]; % near cognate GAA mM^-1s^-1 
kcatKMnc = kcatKMnc/1000; %convert to uM^-1s^-1
errorc = [4; 9; 9; 7; 14; 18]; % error of cognate AAA uM^-1s^-1
errornc = [1; 2; 4; 10; 16; 160]; % error of near cognate GAA mM^-1s^-1 
errornc = errornc/1000; %convert to uM^-1s^-1

kcatKMpepnc = [3.9e-4; 2.7e-3; 9.86e-3;3.67e-2; 2.5e-1];  %don't have Mg = 1.9 , unit is uM^-1s^-1
errorpepnc = [3e-5; 2.3e-4; 1.61e-3; 6.79e-3; 7.18e-2];  % unit is uM^-1s^-1
kcatKMpepc = [60; 117; 147; 167; 180]; %don't have Mg = 1.9

kpepnc = 0.3 ;
%to correspond with qnc 
kcatKMc(2) = []; 
kcatKMnc(2) = [];
errorc(2) = [];
errornc(2) = [];

Jtotal = NaN*ones(81,length(Mgpep));
Atotal = NaN*ones(81,length(Mgpep));
JIS = NaN*ones(81,length(Mgpep));
Q = NaN*ones(81,length(Mgpep));
AIS = NaN*ones(81,length(Mgpep));
Aproof = NaN*ones(81,length(Mgpep));
Jtotalc = NaN*ones(81,length(Mgpep));
Jtotalnc = NaN*ones(81,length(Mgpep));
Qc = NaN*ones(81,length(Mgpep));
Qnc = NaN*ones(81,length(Mgpep));
JISc = NaN*ones(81,length(Mgpep));
JISnc = NaN*ones(81,length(Mgpep));

for i = 1 : length(Mgpep)
    
    %create list of kcatKMc and qnc
    kcatKMc_list = NaN*ones(3,1);
    kcatKMc_list = [kcatKMc(i,1); kcatKMc(i,1) + errorc(i,1); kcatKMc(i,1) - errorc(i,1)];
    kcatKMnc_list = NaN*ones(3,1);
    kcatKMnc_list = [kcatKMnc(i,1); kcatKMnc(i,1) + errornc(i,1); kcatKMnc(i,1) - errornc(i,1)]; 
    kcatKMpepnc_list = NaN*ones(3,1);
    kcatKMpepnc_list = [kcatKMpepnc(i,1); kcatKMpepnc(i,1) + errorpepnc(i,1); kcatKMpepnc(i,1) - errorpepnc(i,1)]; 

    qnc_list = NaN*ones(9,1);
    
    for m = 1 : length(kcatKMnc_list)
        qnc_list(3*(m-1)+1:3*(m-1)+3,1) = (kcatKMnc_list(m,1)./kcatKMpepnc_list-1)*kpepnc;   
    end
    
    [package, dist] = proteincombo(kcatKMc_list, kcatKMnc_list, qnc_list);
    
    Jtotal(:,i) = package{1,1};
    Atotal(:,i) = package{2,1};
    JIS(:,i) =  package{3,1};
    Q(:,i) = package{4,1};
    AIS(:,i) = package{5,1};
    Aproof(:,i) = package{6,1};
    Jtotalc(:,i) = package{7,1};
    Jtotalnc(:,i) = package{8,1};
    Qc(:,i) = package{9,1};
    Qnc(:,i) = package{10,1};
    JISc(:,i) =  package{11,1};
    JISnc(:,i) =  package{12,1}; 
    
    
end

Atotalmedian = NaN*ones(length(Mgpep),1);
Jtotalmedian = NaN*ones(length(Mgpep),1);
lerrx = NaN*ones(length(Mgpep),1);
rerrx = NaN*ones(length(Mgpep),1);
derry = NaN*ones(length(Mgpep),1);
uerry = NaN*ones(length(Mgpep),1);

for i = 1 : length(Mgpep)
    Atotalmedian(i,1) = median(Atotal(:,i));
    Jtotalmedian(i,1) = median(Jtotal(:,i));
    temp = quantile(Atotal(:,i),[0.25 0.75]);
    lerrx(i,1) = Atotalmedian(i,1) - temp(1);rerrx(i,1) = temp(2) - Atotalmedian(i,1);
    temp2 = quantile(Jtotal(:,i),[0.25 0.75]);
    derry(i,1) = Jtotalmedian(i,1) - temp2(1); uerry(i,1) = temp2(2) - Jtotalmedian(i,1);
end

h = boxplot(Jtotal, Atotalmedian );
ax = gca;
ax.XTickLabel = Mgpep;
xlabel('Mg^{2+} (mM)');
ylabel('Normalized product flux')

boxplot(Atotal, Jtotalmedian, 'Orientation', 'horizontal');

errorbarxy(Atotalmedian, Jtotalmedian, lerrx, rerrx, derry, uerry) 
ax = gca;
ax.Box = 'off';
ax.XTickLabel
xlabel('A_{total}')
ylabel('Normalized product flux')

ax2 = axes('Position',ax.Position,'XAxisLocation','top', 'YAxisLocation','left','Color','none','YTick',[],'YTickLabel',[])
data2_tick = Mgpep; 
ax2.XLim = ax.XLim;
ax2.XTick = sort(Atotalmedian,'ascend');
ax2.XTickLabel = sort(Mgpep,'descend');
xlabel('Free Mg^{2+} (mM)')

figure
Mgpeperror = double.empty(length(Mgpep),0);
errorbarxy(Mgpep, Jtotalmedian,Mgpeperror,Mgpeperror,derry, uerry,{'bs-','','b'})
ax = gca;
set(ax,'Box','off','YColor','b');
ylabel('Normalized product flux');
xlabel('Free Mg^{2+} (mM)')

ax2 = axes('Position',ax.Position,'XAxisLocation','bottom','Color','none','XTick',[], 'YAxisLocation','right')
hold(ax2)
errorbarxy(Mgpep, Atotalmedian,Mgpeperror,Mgpeperror,lerrx, rerrx,{'ro-','','r'})
set(ax2,'YScale','log','YColor','r');
ylabel('Total accuracy');
