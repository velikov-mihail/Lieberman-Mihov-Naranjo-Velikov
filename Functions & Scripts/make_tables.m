%% Start a diary (will have all the table results)
clear
clc

diary Results\LNMVTablesOutput.txt

fprintf('\n\n\n\nTable printing started @ %s\n\n\n',char(datetime('now')));

%% Table 2: Paydex distribution

fprintf('\n\n\nTable 2 output:\n\n\n');

clear
clc

load paydex
load dates

% Find the range
s = find(dates==200601);
e = find(dates==201912);

% Paydex distribution
% Panel A - full sample
temp = [nan mean(mean(paydex(s:e,:), 2, 'omitnan'), 'omitnan') ...
            mean(std(paydex(s:e,:), [], 2, 'omitnan'), 'omitnan') ...
            mean(prctile(paydex(s:e,:), [5 25 50 75 95], 2), 'omitnan')];

% Panel B - by year
for i=2006:2019
    ind = find(floor(dates/100)==i);
    temp = [temp; [i mean(mean(paydex(ind,:), 'omitnan'), 'omitnan') ...
                     mean(std(paydex(ind,:), [], 2, 'omitnan'), 'omitnan') ...
                     mean(prctile(paydex(ind,:), [5 25 50 75 95], 2), 'omitnan')]];
end

% Get the row names
h = cellstr(num2str(temp(:,1)));
h(1) = {'PAYDEX'};
temp = temp(:,2:end);

% Print
mat2Tex(temp, temp, h, 2);

%% Table 3: Z_paydex distribution

fprintf('\n\n\nTable 3 output:\n\n\n');

clear
clc

load z_paydex
load dates

s = find(dates==200601);
e = find(dates==201912);

% Z_paydex distribution
% Panel A - full sample
temp = [nan mean(mean(z_paydex(s:e,:), 2, 'omitnan'), 'omitnan') ...
            mean(std(z_paydex(s:e,:), [], 2, 'omitnan'), 'omitnan') ...
            mean(prctile(z_paydex(s:e,:), [5 25 50 75 95], 2), 'omitnan')];

% Panel B - by year
for i=2006:2019
    ind = find(floor(dates/100)==i);
    temp = [temp; [i mean(mean(z_paydex(ind,:), 'omitnan'), 'omitnan') ...
                     mean(std(z_paydex(ind,:), [], 2, 'omitnan'), 'omitnan') ...
                     mean(prctile(z_paydex(ind,:), [5 25 50 75 95], 2), 'omitnan')]];
end

% Get the row names
h = cellstr(num2str(temp(:,1)));
h(1) = {'$Z_{\text{paydex}}$'};
temp = temp(:,2:end);

% Print
mat2Tex(temp, temp, h, 2);

%% Table 5: LPC distribution table

fprintf('\n\n\nTable 5 output:\n\n\n');

clear
clc

load lpc
load dates

s = find(dates==200601);
e = find(dates==201912);

% LPC distribution
% Panel A - full sample
temp = [nan mean(mean(lpc(s:e,:), 2, 'omitnan'), 'omitnan') ...
            mean(std(lpc(s:e,:), [], 2, 'omitnan'), 'omitnan') ...
            mean(prctile(lpc(s:e,:), [5 25 50 75 95], 2), 'omitnan')];

% Panel B - by year
for i=2006:2019
    ind = find(floor(dates/100)==i);
    temp = [temp; [i mean(mean(lpc(ind,:), 'omitnan'), 'omitnan') ...
                     mean(std(lpc(ind,:), [], 2, 'omitnan'), 'omitnan') ...
                     mean(prctile(lpc(ind,:), [5 25 50 75 95], 2), 'omitnan')]];
end

% Get the row names
h = cellstr(num2str(temp(:,1)));
h(1) = {'LPC'};
temp = temp(:,2:end);

% Print
mat2Tex(temp, temp, h, 2);

%% Tables 7, 8, and 9: Z_paydex, Z_slow, and Z_current performance


clear
clc

load z_paydex
load z_d_slow
load z_d_cur

fprintf('\n\n\nTable 7 output:\n\n\n');
printFullSamplePerformanceTable(z_paydex, [200512 201912], 5, 1);

fprintf('\n\n\nTable 8 output:\n\n\n');
printFullSamplePerformanceTable(z_d_slow, [200512 201912], 5, 1);

fprintf('\n\n\nTable 9 output:\n\n\n');
printFullSamplePerformanceTable(z_d_cur, [200512 201912], 5, 1);


%% Table 10: LPC sort

fprintf('\n\n\nTable 10, Panel A+B output:\n\n\n');

clear
clc

load lpc
load dates
load ret
load me
load ff
load NYSE

printFullSamplePerformanceTable(lpc, [200512 201912], 5, 1);

% Sort based on LPC
ind = makeUnivSortInd(lpc, 5, NYSE);
s = find(dates==200601);
e = find(dates==201912);

fprintf('\n\n\nTable 10, Panel A+C output:\n\n\n');


% Controlling for IMC
load IMC
res1 = runUnivSort(ret, ind, dates, me, 'timePeriod', [200512 201912], ...
                                        'plotFigure', 0, ...
                                        'printResults', 0, ...
                                        'factorModel', [mkt(s-1:e) IMC(s-1:e)]);
a = [res1.alpha'; res1.factorLoadings(1).b'; res1.factorLoadings(2).b';];
tA = [res1.talpha'; res1.factorLoadings(1).t'; res1.factorLoadings(2).t';];
heads = [{'$\alpha^{IMC}$'},{'$\beta_{\text{MKT}}$'},{'$\beta_{\text{IMC}}$'}];
mat2Tex(a, tA, heads, 2);




fprintf('\n\n\nTable 10, Panel A+D output:\n\n\n');

% Controlling for CPR
fileURL = 'https://fotig.com/files/research/CPRData_2020.csv';
opts = detectImportOptions(fileURL);
data = readtable(fileURL, opts);

ghs = nan(size(mkt));

[~, ia, ib] = intersect(dates, data.Date);
ghs(ia) = data.R_SSpread(ib)/100;
res1 = runUnivSort(ret, ind, dates, me, 'timePeriod', [200512 201912], ...
                                        'plotFigure', 0, ...
                                        'printResults', 0, ...
                                        'factorModel', [mkt(s-1:e) ghs(s-1:e)]);
a = [res1.alpha'; res1.factorLoadings(1).b'; res1.factorLoadings(2).b';];
tA = [res1.talpha'; res1.factorLoadings(1).t'; res1.factorLoadings(2).t';];
heads = [{'$\alpha^{CR}$'},{'$\beta_{\text{MKT}}$'},{'$\beta_{\text{GHS}}$'}];
mat2Tex(a, tA, heads, 2);

fprintf('\n\n\nTable 10, Panel A+E output:\n\n\n');
% Controlling for industry adjusted CPR
fileURL = 'https://fotig.com/files/research/CPRData_2020.csv';
opts = detectImportOptions(fileURL);
data = readtable(fileURL, opts);

ghs = nan(size(mkt));

[~, ia, ib] = intersect(dates, data.Date);
ghs(ia) = data.R_SSpread_IndustryAdjusted(ib)/100;
res1 = runUnivSort(ret, ind, dates, me, 'timePeriod', [200512 201912], ...
                                        'plotFigure', 0, ...
                                        'printResults', 0, ...
                                        'factorModel', [mkt(s-1:e) ghs(s-1:e)]);
a = [res1.alpha'; res1.factorLoadings(1).b'; res1.factorLoadings(2).b';];
tA = [res1.talpha'; res1.factorLoadings(1).t'; res1.factorLoadings(2).t';];
heads = [{'$\alpha^{adjCR}$'},{'$\beta_{\text{MKT}}$'},{'$\beta_{\text{adjGHS}}$'}];
mat2Tex(a, tA, heads, 2);


%% Table 11

fprintf('\n\n\nTable 11, Fama-MacBeth:\n\n\n');

clear
clc

% Load data
load ret
load dates
load me
load bm
load gp
load R
load OC
load XINT
load AT
load lpc
load lpc_0_50
load lpc_80_100
load iaLPC
load FinFirms
load upStreamness

% Operating leverage
opLev = FillMonths(OC./AT);
opLev(FinFirms==1) = nan;

% Operating leverage with interest expenses
OCstar = nanmatsum(OC, XINT);
opLevStar = FillMonths(OCstar./AT);
opLevStar(FinFirms==1) = nan;

% Gross profitability
gp = GP./AT;
gp(FinFirms==1) = nan;

% Asset growth
ag = AT./lag(AT,12,nan)-1;
ag(FinFirms==1) = nan;

% Kick off nonpositive B/M firms
bm(bm<=0)=nan;

% Cash conversion cycle
load APQ
load INVTQ
load RECTQ
load COGSQ
load SALEQ

avgInv = (INVTQ+lag(INVTQ,3,nan))/2;
avgRec = (RECTQ+lag(RECTQ,3,nan))/2;
avgPay = (APQ+lag(APQ,3,nan))/2;
ccc = 365*(avgInv./(COGSQ+lag(COGSQ,3,nan)) + avgRec./(SALEQ+lag(SALEQ,3,nan)) + avgPay./(COGSQ+lag(COGSQ,3,nan)));
ccc(FinFirms==1) = nan;
ccc = ccc/100;

% Accruals
accruals = getAnomalySignals('novyMarxVelikovAnomalies.csv', 'permno', 'dates', ...
                                                'anomalyNames', {'accruals'});

% Run Fama-MacBeth regressions
res1 = runFamaMacBeth(100*ret, [log(me) log(bm) gp ag R ret], dates, 'timePeriod', [200512 201912], 'printResults', 0);
res2 = runFamaMacBeth(100*ret, [lpc log(me) log(bm) gp ag R ret], dates, 'timePeriod', [200512 201912], 'printResults', 0);
res3 = runFamaMacBeth(100*ret, [lpc_0_50 log(me) log(bm) gp ag R ret], dates, 'timePeriod', [200512 201912], 'printResults', 0);
res4 = runFamaMacBeth(100*ret, [lpc_80_100 log(me) log(bm) gp ag R ret], dates, 'timePeriod', [200512 201912], 'printResults', 0);
res5 = runFamaMacBeth(100*ret, [opLev log(me) log(bm) gp ag R ret], dates, 'timePeriod', [200512 201912], 'printResults', 0);
res6 = runFamaMacBeth(100*ret, [lpc opLev log(me) log(bm) gp ag R ret], dates, 'timePeriod', [200512 201912], 'printResults', 0);
res7 = runFamaMacBeth(100*ret, [opLevStar log(me) log(bm) gp ag R ret], dates, 'timePeriod', [200512 201912], 'printResults', 0);
res8 = runFamaMacBeth(100*ret, [lpc opLevStar log(me) log(bm) gp ag R ret], dates, 'timePeriod', [200512 201912], 'printResults', 0);
res9 = runFamaMacBeth(100*ret, [upStreamness log(me) log(bm) gp ag R ret], dates, 'timePeriod', [200512 201912], 'printResults', 0);
res10 = runFamaMacBeth(100*ret, [lpc upStreamness log(me) log(bm) gp ag R ret], dates, 'timePeriod', [200512 201912], 'printResults', 0);
res11 = runFamaMacBeth(100*ret, [ccc log(me) log(bm) gp ag R ret], dates, 'timePeriod', [200512 201912], 'printResults', 0);
res12 = runFamaMacBeth(100*ret, [lpc ccc log(me) log(bm) gp ag R ret], dates, 'timePeriod', [200512 201912], 'printResults', 0);
res13 = runFamaMacBeth(100*ret, [accruals log(me) log(bm) gp ag R ret], dates, 'timePeriod', [200512 201912], 'printResults', 0);
res14 = runFamaMacBeth(100*ret, [lpc accruals log(me) log(bm) gp ag R ret], dates, 'timePeriod', [200512 201912], 'printResults', 0);
res15 = runFamaMacBeth(100*ret, [iaLPC log(me) log(bm) gp ag R ret], dates, 'timePeriod', [200512 201912], 'printResults', 0);
res16 = runFamaMacBeth(100*ret, [lpc opLev upStreamness ccc accruals log(me) log(bm) gp ag R ret], dates, 'timePeriod', [200512 201912], 'printResults', 0);
res17 = runFamaMacBeth(100*ret, [iaLPC opLev upStreamness ccc accruals log(me) log(bm) gp ag R ret], dates, 'timePeriod', [200512 201912], 'printResults', 0);

a = [res1.bhat(1) nan(1,9) res1.bhat(2:end) 100*res1.mean_R2 sum(isfinite(res1.beta(:,1))); 
     res2.bhat(1:2) nan(1,8) res2.bhat(3:end)  100*res2.mean_R2 sum(isfinite(res2.beta(:,1)));
     res3.bhat(1) nan(1,2) res3.bhat(2) nan(1,6) res3.bhat(3:end)  100*res3.mean_R2 sum(isfinite(res3.beta(:,1)));
     res4.bhat(1) nan(1,3) res4.bhat(2) nan(1,5) res4.bhat(3:end) 100*res4.mean_R2 sum(isfinite(res4.beta(:,1)));
     res5.bhat(1) nan(1,4) res5.bhat(2) nan(1,4) res5.bhat(3:end)  100*res5.mean_R2 sum(isfinite(res5.beta(:,1)));
     res6.bhat(1:2) nan(1,3) res6.bhat(3) nan(1,4) res6.bhat(4:end) 100*res6.mean_R2 sum(isfinite(res6.beta(:,1)));
     res7.bhat(1)  nan(1,5)  res7.bhat(2) nan(1,3) res7.bhat(3:end) 100*res7.mean_R2 sum(isfinite(res7.beta(:,1)));
     res8.bhat(1:2) nan(1,4) res8.bhat(3) nan(1,3) res8.bhat(4:end) 100*res8.mean_R2 sum(isfinite(res8.beta(:,1)));
     res9.bhat(1) nan(1,6)  res9.bhat(2) nan(1,2) res9.bhat(3:end) 100*res9.mean_R2 sum(isfinite(res9.beta(:,1)));
     res10.bhat(1:2) nan(1,5) res10.bhat(3) nan(1,2) res10.bhat(4:end) 100*res10.mean_R2 sum(isfinite(res10.beta(:,1)));
     res11.bhat(1) nan(1,7) res11.bhat(2) nan res11.bhat(3:end) 100*res11.mean_R2 sum(isfinite(res11.beta(:,1)));
     res12.bhat(1:2) nan(1,6) res12.bhat(3) nan res12.bhat(4:end) 100*res12.mean_R2 sum(isfinite(res12.beta(:,1)));
     res13.bhat(1) nan(1,8) res13.bhat(2:end) 100*res13.mean_R2 sum(isfinite(res13.beta(:,1)));
     res14.bhat(1:2) nan(1,7) res14.bhat(3:end) 100*res14.mean_R2 sum(isfinite(res14.beta(:,1)));
     res15.bhat(1) nan res15.bhat(2) nan(1,7) res15.bhat(3:end)  100*res15.mean_R2 sum(isfinite(res15.beta(:,1)));
     res16.bhat(1:2) nan(1,3) res16.bhat(3) nan res16.bhat(4:end)  100*res16.mean_R2 sum(isfinite(res16.beta(:,1)));
     res17.bhat(1) nan res17.bhat(2) nan(1,2) res17.bhat(3) nan res17.bhat(4:end)  100*res17.mean_R2 sum(isfinite(res17.beta(:,1)));
]';
  
% Prepare t-statistics matrix
tA =   [res1.t(1) nan(1,9) res1.t(2:end) nan nan; 
     res2.t(1:2) nan(1,8) res2.t(3:end) nan nan;
     res3.t(1) nan(1,2) res3.t(2) nan(1,6) res3.t(3:end) nan nan;
     res4.t(1) nan(1,3) res4.t(2) nan(1,5) res4.t(3:end) nan nan;
     res5.t(1) nan(1,4) res5.t(2) nan(1,4) res5.t(3:end) nan nan;
     res6.t(1:2) nan(1,3) res6.t(3) nan(1,4) res6.t(4:end) nan nan;
     res7.t(1)  nan(1,5)  res7.t(2) nan(1,3) res7.t(3:end) nan nan;
     res8.t(1:2) nan(1,4) res8.t(3) nan(1,3) res8.t(4:end) nan nan;
     res9.t(1) nan(1,6)  res9.t(2) nan(1,2) res9.t(3:end) nan nan;
     res10.t(1:2) nan(1,5) res10.t(3) nan(1,2) res10.t(4:end) nan nan;
     res11.t(1) nan(1,7) res11.t(2) nan res11.t(3:end) nan nan;
     res12.t(1:2) nan(1,6) res12.t(3) nan res12.t(4:end) nan nan;
     res13.t(1) nan(1,8) res13.t(2:end) nan nan;
     res14.t(1:2) nan(1,7) res14.t(3:end) nan nan;
     res15.t(1) nan res15.t(2) nan(1,7) res15.t(3:end) nan nan;
     res16.t(1:2) nan(1,3) res16.t(3) nan res16.t(4:end) nan nan;
     res17.t(1) nan res17.t(2) nan(1,2) res17.t(3) nan res17.t(4:end) nan nan;
]';

% Define header
h = {'Const','LPC','LPC$^{adj.}$','LPC$_{0-50}$','LPC$_{80-100}$','Operating Leverage','Operating Leverage$^+$','Upstreamness','CCC','Accruals', ...
     'log(ME)','log(B/M)','GP/A','AG','$r_{12,1}$','$r_{1,0}$','$R^2 (\%)$','n'};

% Generate LaTeX table
mat2Tex(a, tA, h, 2);

% Generate separate row for number of observations
mat2Tex(a(end,:), a(end,:), {'$n$'}, 0);



%% Table 12: GMM

fprintf('\n\n\nTable 11, GMM:\n\n\n');

clear
clc

load dates
load lpc
load dates
load nyse
load ret
load me
load ff
load factors

% Sort on LPC into deciles
ind = makeUnivSortInd(lpc, 10);
res = runUnivSort(ret, ind, dates, me, 'plotFigure', 0, ...
                                       'printResults', 0);

% Extract the test asset returns
testAssets = res.pret(:,1:end-1);
s = find(dates==200601);
e = find(dates==201912);
testAssets(1:s-1, :) = nan;
testAssets(e+1, :) = nan;

% Store the number of factors
nFactors = length(factors);

clear res

% Loop through the factors & estimate GMM
for i=1:(nFactors+1)
    if i<=nFactors
        x = [mkt factors(i).factor];   
    else
        x = mkt;
    end
    ind = isfinite(sum([x testAssets],2));

    [tempRes] = runGMM(x(ind,:), rf(ind), testAssets(ind,:), 1);   
    res(i,1).res = tempRes;
end

% Print the results
tablePtfs = [1 5 10];

labels = [{'MKT'} factors.label];
labels = regexprep(labels,'_','\\_');

panelA = [100*res(end).res.alphas([tablePtfs end]); res(end).res.talphas(end);];
panelB = [nan(5,1)];

fprintf(' & MKTRF');
for i=1:nFactors
    fprintf(' & %s', char(factors(i).label));
    panelA = [panelA [100*res(i).res.alphas([tablePtfs end]); res(i).res.talphas(end);]];
    panelB = [panelB [100*res(i).res.betas([tablePtfs end]); res(i).res.tbetas(end)]];
end
fprintf('\\\\[2pt]');

h = {'Low LPC (L)','Medium LPC','High LPC (H)'};
mat2Tex(panelA(1:3,:), panelA(1:3,:), h, 2);

h = {'Spread (H-L)'};
mat2Tex(panelA(4,:), panelA(5,:), h, 2);


h = {'Low LPC (L)','Medium LPC','High LPC (H)'};
mat2Tex(panelB(1:3,:), panelB(1:3,:), h, 2);

h = {'Spread (H-L)'};
mat2Tex(panelB(4,:), panelB(5,:), h, 2);

% Panel C
% 5 industries
fileURL = 'https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/5_Industry_Portfolios_CSV.zip';
ffFileName = unzip(fileURL, 'Data/');

% Read in the 3 factors
opts = detectImportOptions(char(ffFileName));
FFptfs = readtable(char(ffFileName), opts);

% Clean up the file - it also includes annual returns for the factors;
e = find(isnan(FFptfs.Var1), 1, 'first');
FFptfs(e:end,:) = [];


[~, ia, ib] = intersect(dates, FFptfs.Var1);
ffPtfs = nan(length(dates), size(FFptfs, 2)-1);
ffPtfs(ia, :) = table2array(FFptfs(ib, 2:end))/100;
testAssets = [testAssets ffPtfs];


% 6 size - btm portfolios
fileURL = 'https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/6_Portfolios_2x3_CSV.zip';
ffFileName = unzip(fileURL, 'Data/');

% Read in the 3 factors
opts = detectImportOptions(char(ffFileName));
FFptfs = readtable(char(ffFileName), opts);

% Clean up the file - it also includes annual returns for the factors;
e = find(isnan(FFptfs.Var1), 1, 'first');
FFptfs(e:end,:) = [];


[~, ia, ib] = intersect(dates, FFptfs.Var1);
ffPtfs = nan(length(dates), size(FFptfs, 2)-1);
ffPtfs(ia, :) = table2array(FFptfs(ib, 2:end))/100;
testAssets = [testAssets ffPtfs];

% 6 size - prof portfolios
fileURL = 'https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/6_Portfolios_ME_OP_2x3_CSV.zip';
ffFileName = unzip(fileURL, 'Data/');

% Read in the 3 factors
opts = detectImportOptions(char(ffFileName));
FFptfs = readtable(char(ffFileName), opts);

% Clean up the file - it also includes annual returns for the factors;
e = find(isnan(FFptfs.Var1), 1, 'first');
FFptfs(e:end,:) = [];


[~, ia, ib] = intersect(dates, FFptfs.Var1);
ffPtfs = nan(length(dates), size(FFptfs, 2)-1);
ffPtfs(ia, :) = table2array(FFptfs(ib, 2:end))/100;
testAssets = [testAssets ffPtfs];



% 6 size - inv portfolios
fileURL = 'https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/6_Portfolios_ME_INV_2x3_CSV.zip';
ffFileName = unzip(fileURL, 'Data/');

% Read in the 3 factors
opts = detectImportOptions(char(ffFileName));
FFptfs = readtable(char(ffFileName), opts);

% Clean up the file - it also includes annual returns for the factors;
e = find(isnan(FFptfs.Var1), 1, 'first');
FFptfs(e:end,:) = [];


[~, ia, ib] = intersect(dates, FFptfs.Var1);
ffPtfs = nan(length(dates), size(FFptfs, 2)-1);
ffPtfs(ia, :) = table2array(FFptfs(ib, 2:end))/100;
testAssets = [testAssets ffPtfs];


% Loop through the factors & estimate GMM
for i=1:(nFactors+1)
    if i<=nFactors
        x = [mkt factors(i).factor];   
    else
        x = mkt;
    end
    ind = isfinite(sum([x testAssets],2));

    [tempRes] = runGMM(x(ind,:), rf(ind), testAssets(ind,:), 1);   
    res2(i,1).res = tempRes;
end



h = {'$b_{MKTRF}','$b_{GSCPI}$','MAE'};
a = [res(end).res.theta(1) res(1).res.theta(1) nan(1,3) res2(end).res.theta(1) res2(1).res.theta(1) ;
     nan res(1).res.theta(2) nan(1,3) nan res2(1).res.theta(2);
     res(end).res.MAE res(1).res.MAE nan(1,3) res2(end).res.MAE res2(1).res.MAE];
tA = [res(end).res.t(1) res(1).res.t(1) nan(1,3) res2(end).res.t(1) res2(1).res.t(1);
     nan res(1).res.t(2) nan(1,3) nan res2(1).res.t(2);
     nan nan nan(1,3) nan nan];
mat2Tex(a, tA, h);




diary('off');
