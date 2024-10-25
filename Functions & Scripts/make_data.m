%% Make variables from CSAD

clear
clc

load permno
load dates
load ret

% Read the input CSAD (Credit Score Archive Database) data from D&B
opts = detectImportOptions('CRSP_DNB_CSAD.csv');
data = readtable('CRSP_DNB_CSAD.csv',opts);

% Store a few constants
[nMonths, nStocks] = size(ret);
nObs = height(data);

% Initialize matrices for the main variables of interest
paydex = nan(nMonths, nStocks);
d_slow = nan(nMonths, nStocks);
d_cur = nan(nMonths, nStocks);

% Loop through the CSAD data and assign the values to the matrices
for i=1:nObs
    c = find(permno==data.PERMNO(i));
    r = find(dates==data.Actual_Archive_Date(i));
    if isfinite(r+c)
        paydex(r,c) = data.pydx_1(i);        
        d_slow(r,c) = data.d_slow(i);
        d_cur(r,c) = data.d_cur(i);
    end
end    

% Store the variables
save Data\paydex paydex
save Data\d_cur d_cur
save Data\d_slow d_slow

%% Make standardized variables 

clear
clc

load ret
load nyse
load dates
load permno
load me
load paydex
load d_slow
load d_cur

% We'll start in 2005
s = find(dates==200501);

% Store a few constants
[nMonths, nStocks] = size(ret);

% Initialize mean & std of paydex matrices
mean_paydex = nan(nMonths, nStocks);
std_paydex = nan(nMonths, nStocks);

% Loop through the months
for i=s:nMonths
    
    % We'll require at least 3 observations
    ind = sum(isfinite(paydex(i-11:i,:)),1) >= 3;

    % Calculate the average paydex
    mean_paydex(i,ind) = mean(paydex(i-11:i,ind), 'omitnan');

    % Std will be 1st x-sectional percentile if zoro
    denom = std(paydex(i-11:i,ind), [], 1, 'omitnan');
    denom(denom==0) = prctile(denom(denom~=0), 1);
    std_paydex(i,ind) = denom;
end

% Calculate z-paydex
z_paydex = (paydex-lag(mean_paydex,1,nan))./lag(std_paydex,1,nan);

% Initialize the d_slow mean & std variables
mean_d_slow = nan(nMonths, nStocks);
std_d_slow = nan(nMonths, nStocks);

% Loop through the months
for i=s:nMonths

    % We'll require at least 3 observations
    ind = sum(isfinite(d_slow(i-11:i,:)), 1) >= 3;

    % Calculate the average d_slow
    mean_d_slow(i,ind) = mean(d_slow(i-11:i,ind), 'omitnan');

    % Std will be 1st x-sectional percentile if zoro
    denom = std(d_slow(i-11:i,ind), [], 1, 'omitnan');    
    denom(denom==0) = prctile(denom(denom~=0), 1);
    std_d_slow(i,ind) = denom;
end

% Calculate z-dollars-slow
z_d_slow = (d_slow-lag(mean_d_slow,1,nan))./lag(std_d_slow,1,nan);

% Initialize the d_cur mean & std variables
mean_d_cur = nan(nMonths, nStocks);
std_d_cur = nan(nMonths, nStocks);

% Loop through the months
for i=s:nMonths   
    % We'll require at least 3 observations
    ind = sum(isfinite(d_cur(i-11:i,:)),1)>=3;

    % Calculate the average d_cur
    mean_d_cur(i,ind)=mean(d_cur(i-11:i,ind), 'omitnan');

    % Std will be 1st x-sectional percentile if zoro
    denom = std(d_cur(i-11:i,ind), [] ,1, 'omitnan');
    denom(denom==0) = prctile(denom(denom~=0), 1);
    std_d_cur(i,ind) = denom;
end

% Calculate z-dollars-current
z_d_cur = (d_cur-lag(mean_d_cur,1,nan))./lag(std_d_cur,1,nan);

% Store the variables
save Data\z_paydex z_paydex
save Data\z_d_slow z_d_slow
save Data\z_d_cur z_d_cur

%% Make LPC 

clear
clc

load ret
load nyse
load dates
load permno
load me
load paydex
load ff

% We'll start in 2005
s = find(dates==200501);

% Store a few constants
[nMonths, nStocks] = size(ret);

% Initialize the LPC structure
lpcStruct = struct;

% Loop through the rolling windows
for nroll = 12:12:120
    % initialize the LPC for this rolling window
    lpc = nan(nMonths, nStocks);
    
    % Loop through the months
    for i=s:nMonths 

        % Calculate LPC
        lpc(i,:) = sum(paydex(i-nroll+1:i, :) >= 70 & ...
                       paydex(i-nroll+1:i, :) < 80, 1) ./ ...
                   sum(isfinite(paydex(i-nroll+1:i, :)));
    end   
    
    % Set observations for which we don't have paydex to NaN
    lpc(isnan(paydex)) = nan;
    
    % Store the LPC & rolling window
    count = nroll/12;
    lpcStruct(count,1).lpc = lpc;
    lpcStruct(count,1).nroll = nroll;
end

% Test the LPC
lpc = lpcStruct(5).lpc;
ind = makeUnivSortInd(lpc, 5, NYSE);
res1 = runUnivSort(ret, ind, dates, me, 'timePeriod', 200512, ...
                                        'factorModel', 6, ...
                                        'plotFigure', 0);

% Save both the 5-year LPC and the structure
save Data\lpc lpc
save Data\lpcStruct lpcStruct -v7.3

%% Download a few more COMPUSTAT variables

clear
clc

% getCOMPUSTATAdditionalData(usernameUI(), passwordUI(), {'APQ', 'RECTRQ'}, 'quarterly');

%% Make WHited-Wu

clear
clc

% Load the data
load ATQ
load SALEQ
load IBQ
load DLTTQ
load SIC
load DVC
load DPQ
load DVP
load dates
load permno
load FQTR
load RDQ
load DATAFQTR
load LPC

% Drop the repeating RDQ and FQTR observations
indToDrop = (RDQ - lag(RDQ, 1, nan)) == 0;
RDQ(indToDrop)=nan;
FQTR(indToDrop)=nan;
DATAFQTR(indToDrop)=nan;

% We divide/take logs of these two
ATQ(ATQ<=0) = nan;
SALEQ(SALEQ<=0) = nan;

% We are using annual dividends as there is no quarterly DVCQ
dividends = FillMonths(DVC + DVP);
posDiv = 1 * (dividends > 0);
posDiv(~isfinite(ATQ)) = nan;

% Prepare some variables
cashFlow = (IBQ + DPQ)./ATQ;
totalDebt = DLTTQ./ATQ;
logAssets = log(ATQ);

SIC(SIC==0) = nan;
SIC = floor(SIC/10);
RDQ = floor(RDQ/100);

% Reshape & create a table
nStocks = size(ATQ, 2);
nMonths = size(ATQ, 1);
nObs = nStocks * nMonths;
rptdPermno = repmat(permno', nMonths, 1);
data = [reshape(FQTR,  nObs, 1) ...
        reshape(rptdPermno, nObs, 1) ...
        reshape(RDQ,        nObs, 1) ...
        reshape(cashFlow,   nObs, 1) ... 
        reshape(posDiv,     nObs, 1) ...
        reshape(totalDebt,  nObs, 1) ...
        reshape(logAssets,  nObs, 1) ...
        reshape(SIC,        nObs, 1) ...
        reshape(SALEQ,      nObs, 1)];
indToDrop = isnan(data(:,1));
data(indToDrop, :) = [];

% Make the table now
data = array2table(data);
data.Properties.VariableNames = {'FQTR','permno','RDQ','cashFlow','posDiv', ...
                                 'totalDebt','logAssets','SIC','SALEQ'};

% Create the YYYYQQ variable
data.FQTR = datetime(data.FQTR, 'ConvertFrom', 'YYYYmmDD');
data.month = month(data.FQTR);
data.qtr = floor((data.month-1)/3) + 1;
data.yyyyqq = year(data.FQTR)*100 + data.qtr;

% Calculate industry sales for the industry sales growth variable
dataInd = varfun(@(x) sum(x, 'omitnan'), data(:, {'SIC','yyyyqq','SALEQ'}), ...
                                    'GroupingVariables', {'SIC','yyyyqq'}, ...
                                    'InputVariables',{'SALEQ'});
dataInd.GroupCount = [];
dataInd.Properties.VariableNames={'SIC','yyyyqq','indSales'};

% Merge the industry data back to the stock-level data
data = outerjoin(data, dataInd, 'Type', 'Left', ...
                                'MergeKeys', 1);

% Prepare the variables for which we want to take lags
dataToMerge = data(:, {'permno','SALEQ','indSales','FQTR'});
dataToMerge.Properties.VariableNames = {'permno','lagSALEQ','lagIndSale','lagFQTR'};

dataLag = outerjoin(data, dataToMerge, 'Keys', 'permno');
indToDrop = calmonths(between(dataLag.FQTR, dataLag.lagFQTR)) ~= -3;
dataLag(indToDrop, :) = [];
dataLag = dataLag(:,{'permno_data','FQTR','lagFQTR','lagSALEQ','lagIndSale'});
dataLag.Properties.VariableNames={'permno','FQTR','lagFQTR','lagSALEQ','lagIndSale'};

% Attach the lags
data = outerjoin(data, dataLag, 'Type', 'Left', ...
                                'MergeKeys', 1, ...
                                'RightVariables', {'lagSALEQ','lagIndSale'});

% Create the sales and industry sales growth variables
data.salesGrowth = data.SALEQ./data.lagSALEQ - 1;
data.indSalesGrowth = data.indSales./data.lagIndSale - 1;

data = data(:,{'permno','RDQ','cashFlow','posDiv','totalDebt', ...
                'logAssets','indSalesGrowth','salesGrowth'});

data.WW = -0.091 * data.cashFlow + ... 
          -0.062 * data.posDiv + ... 
           0.021 * data.totalDebt + ... 
          -0.044 * data.logAssets + ... 
           0.102 * data.indSalesGrowth + ... 
          -0.035 * data.salesGrowth;

data = sortrows(data, {'permno','RDQ'}, 'ascend');
data = varfun(@(x) x(end,:), data, 'GroupingVariables', {'permno','RDQ'});
data.GroupCount = [];    
data.Properties.VariableNames = regexprep(data.Properties.VariableNames, 'Fun_','');

% Load a few variables
load permno
load dates
load ret

% Store a few constants
nStocks = length(permno);
nMonths = length(dates);
nObs = nStocks * nMonths;

% Create the linking table with CRSP
rptdDates = repmat(dates, 1, nStocks);
rptdPermno = repmat(permno', nMonths, 1);
crspMatLink = [reshape(rptdPermno, nObs, 1) ...
               reshape(rptdDates, nObs, 1)];
crspMatLinkTab = array2table(crspMatLink, 'VariableNames', {'permno', 'dates'});           


mergedTab = outerjoin(crspMatLinkTab, data, 'Type', 'Left', ...
                                            'LeftKeys', {'permno','dates'}, ...
                                            'RightKeys',{'permno','RDQ'}, ...
                                            'RightVariables', 'WW', ...
                                            'MergeKeys', 1);

% Unstack the table and turn it into a matrix
thisVar = unstack(mergedTab, 'WW', 'dates_RDQ');
thisVar.permno = [];
WhitedWu = table2array(thisVar)';

% Fill in the missing months if any
WhitedWu = FillMonths(WhitedWu);
WhitedWu(isnan(ATQ)) = nan;

save Data\WhitedWu WhitedWu

%% Make other variables (for LPC regs)

clear
clc

load ret
load nyse
load dates
load permno
load me
load paydex

% Store a few constants
[nMonths, nStocks] = size(ret);

% Market share
load SALE
load FF49

% Initialize the market share matrix
market_share = nan(nMonths, nStocks);

% Loop through the industries
for i=1:49
    temp = SALE;
    temp(FF49~=i) = nan;
    tempTotalMktShr = repmat(sum(temp, 2, 'omitnan'), 1, nStocks);
    temp_mkt_shr = temp./tempTotalMktShr;
    market_share(FF49==i) = temp_mkt_shr(FF49==i);
end

% Store the market share variable
save Data\market_share market_share


% Calculate price cost margin
load COGS
load XSGA
load SIC
load OIADP
load FinFirms

% Set nonpositive SALE obs to NaN
SALE(SALE<=0) = nan;

% Calculate operating profits
op_profits = nanmatsum(SALE-COGS, -XSGA);
ind = isnan(op_profits) & isfinite(OIADP);
op_profits(ind) = OIADP(ind);

% Start with firm-level price cost margin
fpcm = op_profits./SALE;

% Calculate the industry one at the 2-digit SIC industry

% Start with 2-digit SIC code
sic = floor(SIC/100);

% Get the unique number of SIC codes
usic = unique(sic);
usic(isnan(usic) | usic==0) = [];
nUsic = length(usic);

% Initialize the industry price-cost margin
ipcm = nan(nMonths, nStocks);
s = find(dates==199406);

% Loop through the years
for i=s:12:nMonths
    % Loop through the industries
    for j=1:nUsic
        % Calculate the industry price cost margin as a weighted average of
        % the firm-level one
        ind_row = find(sic(i,:) == usic(j) & ....
                       isfinite(fpcm(i,:)));
        ipcm(i,ind_row) = sum((fpcm(i,ind_row).*SALE(i,ind_row)) / sum(SALE(i,ind_row)));
    end
end

% Caclculate the industry-adjusted PCM and store it
pcm = fpcm - ipcm;
pcm(SIC<2000 | SIC>=6000) = nan;
save Data\pcm pcm


% Get PMF
% Download & unzip the data from the Hoberg-Phillips website
urlLink = 'https://hobergphillips.tuck.dartmouth.edu/idata/FluidityData.zip';
unzip(urlLink, 'Data');

% Read the PMF data in
opts = detectImportOptions('FluidityData.txt');
data = readtable('FluidityData.txt', opts);

% Read the CRSP/COMPUSTAT linking file in
opts = detectImportOptions('crsp_ccmxpf_lnkhist.csv');
crsp_ccmxpf_lnkhist = readtable('crsp_ccmxpf_lnkhist.csv', opts);

% Filter based on the CCM link
% Leave only link types LC or LU
idxToKeep = ismember(crsp_ccmxpf_lnkhist.linktype, {'LC','LU'});
crsp_ccmxpf_lnkhist = crsp_ccmxpf_lnkhist(idxToKeep, :);

% Leave the variables we need
crsp_ccmxpf_lnkhist = crsp_ccmxpf_lnkhist(:,{'lpermno','gvkey','linkdt','linkenddt'});

% Replace the missing linkeddt with the end of the sample
indNatEndDate = isnat(crsp_ccmxpf_lnkhist.linkenddt);
crsp_ccmxpf_lnkhist.linkenddt(indNatEndDate) = datetime(floor(dates(end)/100), 12, 31);
crsp_ccmxpf_lnkhist.linkdt = year(crsp_ccmxpf_lnkhist.linkdt);
crsp_ccmxpf_lnkhist.linkenddt = year(crsp_ccmxpf_lnkhist.linkenddt);

% Merge & clean the COMPUSTAT annual file with CRSP Link history file
data_linked = outerjoin(data, crsp_ccmxpf_lnkhist, 'Type', 'Left', ...
                                                   'Keys', 'gvkey', ...
                                                   'MergeKeys', 1);

% Fiscal period end date must be within link date range 
idxToDrop = data_linked.year > data_linked.linkenddt |  ... 
            data_linked.year < data_linked.linkdt;
data_linked(idxToDrop,:) = [];  

% Drop missing permnos
idxToDrop = isnan(data_linked.lpermno);
data_linked(idxToDrop,:)=[]; 

% Create the date
data_linked.dates = 100*(data_linked.year+1) + 6;
data_linked.Properties.VariableNames{'lpermno'} = 'permno';

% Clean up
data_linked = data_linked(:,{'permno','dates','prodmktfluid'});
data_linked = sortrows(data_linked, {'permno','dates'});
data_linked = varfun(@(x) x(end,:), data_linked,'GroupingVariables',{'permno','dates'});
data_linked.GroupCount = [];

% merge with CRSP
load crsp_link
crsp_link = outerjoin(crsp_link, data_linked, 'Type', 'Left', ...
                                               'MergeKeys', 1);
% Unstack & store
pmf = unstack(crsp_link, 'Fun_prodmktfluid', 'dates');
pmf.permno = [];
pmf = table2array(pmf)';

save Data\pmf pmf


%% Export signals for submission to Assaying Anomalies

clear
clc

load dates
load permno
load lpc
load z_paydex
load z_d_slow

[nMonths, nStocks] = size(lpc);

rptdDates = repmat(dates, 1, nStocks);
rptdPermno = repmat(permno', nMonths, 1);

s = find(dates==200512);
lpc(1:s,:) = nan;
z_paydex(1:s,:) = nan;
z_d_slow(1:s,:) = nan;

n = nMonths * nStocks;

temp = [reshape(rptdPermno, n, 1) ...
        reshape(rptdDates, n, 1) ...
        reshape(lpc, n, 1) ...
        reshape(z_paydex, n, 1) ...
        reshape(-z_d_slow, n, 1)];
    
% Clean up
temp(sum(isnan(temp(:,3:5)), 2) ==3, :) = [];

% Export
temp = array2table(temp, 'VariableNames', {'permno','dates','lpc','zpaydex','zslow'});
writetable(temp, 'Data\sorting_vars.csv');

data = temp(:, {'permno','dates','lpc'});
writetable(data, 'Data\lpc.csv');

data = temp(:, {'permno','dates','zpaydex'});
writetable(data, 'Data\zpaydex.csv');

data = temp(:, {'permno','dates','zslow'});
writetable(data, 'Data\zslow.csv');


%% Make upstreamness measure 

%% Make NAICS from CRSP's MSEEXCHDATES dataset

clear
clc

load permno
load dates
load ret

opts = detectImportOptions('crsp_mseexchdates.csv');
data = readtable('crsp_mseexchdates.csv', opts);

data = data(:, {'permno','namedt','nameendt','naics'});
data.namedt = 100*year(data.namedt) + month(data.namedt);
data.nameendt = 100*year(data.nameendt) + month(data.nameendt);

data(isnan(data.naics), :) = [];

NAICS = nan(size(ret));

for i = 1:height(data)
    c = find(permno==data.permno(i));
    b = find(dates==data.namedt(i));
    e = find(dates==min(data.nameendt(i), dates(end)));
    NAICS(b:e, c) = data.naics(i);
end

save Data\NAICS NAICS

%% Import the BEA upstreamness measure from R concordance package source

clear
clc

% Import the upstreamness data
opts = detectImportOptions('upstream_us_detailed.csv');
opts.VariableTypes(find(strcmp(opts.VariableNames, 'CODE'))) = {'char'};
data = readtable('upstream_us_detailed.csv', opts);

% Import the 2002 BEA->NAICS mapping
opts = detectImportOptions('bea2002_naics2002.csv');
opts.VariableTypes(find(strcmp(opts.VariableNames, 'BEA2002'))) = {'char'};
link2002 = readtable('bea2002_naics2002.csv', opts);
link2002.Properties.VariableNames = {'Var1','CODE','NAICS_6d','NAICS_5d','NAICS_4d','NAICS_3d','NAICS_2d'};
link2002.BEA_CLASS = repmat({'BEA2002'}, height(link2002), 1);

% Import the 2007 BEA->NAICS mapping
opts = detectImportOptions('bea2007_naics2007.csv');
opts.VariableTypes(find(strcmp(opts.VariableNames, 'BEA2007'))) = {'char'};
link2007 = readtable('bea2007_naics2007.csv', opts);
link2007.Properties.VariableNames = {'Var1','CODE','NAICS_6d','NAICS_5d','NAICS_4d','NAICS_3d','NAICS_2d'};
link2007.BEA_CLASS = repmat({'BEA2007'}, height(link2007), 1);

% Import the 2012 BEA->NAICS mapping
opts = detectImportOptions('bea2012_naics2012.csv');
opts.VariableTypes(find(strcmp(opts.VariableNames, 'BEA2012'))) = {'char'};
link2012 = readtable('bea2012_naics2012.csv', opts);
link2012.Properties.VariableNames = {'Var1','CODE','NAICS_6d','NAICS_5d','NAICS_4d','NAICS_3d','NAICS_2d'};
link2012.BEA_CLASS = repmat({'BEA2012'}, height(link2012), 1);

% Combine & merge the links
link = [link2002; link2007; link2012];
link = sortrows(link, {'CODE','BEA_CLASS'});
link = varfun(@(x) x(end,:), link, 'GroupingVariables', {'CODE','BEA_CLASS'});
link(:,{'GroupCount','Fun_Var1'}) = [];
link.Properties.VariableNames = {'CODE','BEA_CLASS','NAICS_6d','NAICS_5d','NAICS_4d','NAICS_3d','NAICS_2d'};

data = outerjoin(data, link, 'Type', 'Left', ...
                             'Keys', {'CODE','BEA_CLASS'}, ...
                             'RightVariables', {'NAICS_6d','NAICS_5d','NAICS_4d','NAICS_3d','NAICS_2d'}, ...
                             'MergeKeys', 1);

% Now let's work on merging the upstreamness data to our library based on
% NAICS code
load NAICS
load dates

% Create the shorter NAICS codes
n5 = nan(size(NAICS));
n5(NAICS<=99999) = NAICS(NAICS<=99999);
n5(NAICS>99999) = floor(NAICS(NAICS>99999)/10);

n4 = nan(size(NAICS));
n4(NAICS>99999) = floor(NAICS(NAICS>99999)/100);
n4(NAICS>9999 & NAICS<=99999) = floor(n4(NAICS>9999 & NAICS<=99999)/10);
n4(NAICS<=9999) = NAICS(NAICS<=9999);

n3 = nan(size(NAICS));
n3(NAICS>99999) = floor(NAICS(NAICS>99999)/1000);
n3(NAICS>9999 & NAICS<=99999) = floor(NAICS(NAICS>9999 & NAICS<=99999)/100);
n3(NAICS>999 & NAICS<=9999) = floor(NAICS(NAICS>999 & NAICS<=9999)/10);
n3(NAICS<=999) = NAICS(NAICS<=999);

n2 = nan(size(NAICS));
n2(NAICS>99999) = floor(NAICS(NAICS>99999)/10000);
n2(NAICS>9999 & NAICS<=99999) = floor(NAICS(NAICS>9999 & NAICS<=99999)/1000);
n2(NAICS>999 & NAICS<=9999) = floor(NAICS(NAICS>999 & NAICS<=9999)/100);
n2(NAICS>99 & NAICS<=999) = floor(NAICS(NAICS>99 & NAICS<=999)/10);
n2(NAICS<=99) = NAICS(NAICS<=99);

% initialize the upstreamness
upStreamness = nan(size(NAICS));

% Store the BEA release years
yrs = [2002 2007 2012];

% Loop over the years
for r = find(dates==200306):12:length(dates)

    % Find the relevant release
    yr = yrs(find(floor(dates(r)/100)>yrs, 1, 'last'));
    
    % Store the data from that release
    tempData = data(data.YEAR==yr,:);

    % Loop & assign based on 6-digit NAICS match
    for i = 1:height(tempData)           
        ind = NAICS(r,:) == tempData.NAICS_6d(i);
        upStreamness(r, ind) = tempData.GVC_Ui(i);
    end

    % Loop & assign based on 5-digit NAICS match
    for i = 1:height(tempData)          
        ind = isnan(upStreamness(r,:)) & ...
              NAICS(r,:) ~= tempData.NAICS_6d(i) & ...
              n5(r,:) == tempData.NAICS_5d(i);
        upStreamness(r, ind) = tempData.GVC_Ui(i);
    end

    % Loop & assign based on 4-digit NAICS match
    for i = 1:height(tempData)           
        ind = isnan(upStreamness(r,:)) & ...
              NAICS(r,:) ~= tempData.NAICS_6d(i) & ...
              n5(r,:) ~= tempData.NAICS_5d(i) & ...
              n4(r,:) == tempData.NAICS_4d(i);
        upStreamness(r, ind) = tempData.GVC_Ui(i);
    end

    % Loop & assign based on 3-digit NAICS match
    for i = 1:height(tempData)                     
        ind = isnan(upStreamness(r,:)) & ...
              NAICS(r,:) ~= tempData.NAICS_6d(i) & ...
              n5(r,:) ~= tempData.NAICS_5d(i) & ...
              n4(r,:) ~= tempData.NAICS_4d(i) & ...
              n3(r,:) == tempData.NAICS_3d(i);
        upStreamness(r, ind) = tempData.GVC_Ui(i);
    end

    % Loop & assign based on 2-digit NAICS match
    for i = 1:height(tempData)           
        ind = isnan(upStreamness(r,:)) & ...
              NAICS(r,:) ~= tempData.NAICS_6d(i) & ...
              n5(r,:) ~= tempData.NAICS_5d(i) & ...
              n4(r,:) ~= tempData.NAICS_4d(i) & ...
              n3(r,:) ~= tempData.NAICS_3d(i) & ...
              n2(r,:) == tempData.NAICS_2d(i);
        upStreamness(r, ind) = tempData.GVC_Ui(i);
    end
    toc
end

% Fill and store the upstreamness measure
upStreamness = FillMonths(upStreamness);
save Data\upStreamness upStreamness

%% Make factors for GMM

clear
clc

load dates

% Start with constructing a structure with the factors
factors = struct;

% GSCPI index
% Download from https://www.newyorkfed.org/medialibrary/research/interactives/gscpi/downloads/gscpi_data.xlsx
opts = detectImportOptions('gscpi.xlsx');
data = readtable('gscpi.xlsx', opts);

% Fix the date format
data.Var1 = datetime(data.Var1);
data.Var1 = 100*year(data.Var1) + month(data.Var1);

% Assign to a vector
[~, ia, ib] = intersect(dates, data.Var1);
gscpi = nan(size(dates));
gscpi(ia) = data.Var2(ib);

% Take the first difference & assign to structure
d_gscpi =  gscpi-lag(gscpi, 1, nan);

factors(1).factor = d_gscpi;
factors(1).label = 'GSCPI';
clearvars -except factors dates

% He, Kelly, and Manela (2017) - Intermediary Capital Ratio and Risk
% Factors
fileUrl = 'https://zhiguohe.net/wp-content/uploads/2024/07/He_Kelly_Manela_Factors_monthly_240723.csv';
opts = detectImportOptions(fileUrl);
data = readtable(fileUrl, opts);
[~, ia, ib] = intersect(dates,data.yyyymm);
factors(2).factor = nan(size(dates));
factors(2).label = 'INT';
factors(2).factor(ia) =  data.intermediary_capital_risk_factor(ib);
clearvars -except factors dates


% Pastor-Stambaugh factor
fileURL = 'https://finance.wharton.upenn.edu/~stambaug/liq_data_1962_2022.txt';
opts = detectImportOptions(fileURL);
data = readtable(fileURL, opts);
[~, ia, ib] = intersect(dates, data{:,1});
factors(3).factor = nan(size(dates));
factors(3).factor(ia) = data.InnovLiq_eq8_(ib);
factors(3).label = {'LIQ'};
clearvars -except factors dates


% Jurado, Ludvigson, and Ng uncertainty
% getGoogleDriveData('MacroFinanceUncertainty_202306.zip', '1LfVHGC1wBoLNv1xCzAMBUbI4UAgJWEgg','Data\');
% unzip('MacroFinanceUncertainty_202306.zip', 'Data\');
fileName = 'RealUncertaintyToCirculate.xlsx';
opts = detectImportOptions(fileName);
data = readtable(fileName, opts);
data.Date = 100*year(data.Date) + month(data.Date);
[~, ia, ib] = intersect(dates, data.Date);
factors(4).factor = nan(size(dates));
temp = data.h_1(ib);
factors(4).factor(ia) = log(temp) - log(lag(temp,1,nan));
factors(4).label = {'UNC'};
clearvars -except factors dates

% Create the aggregate TC factor from Grigoris, Hu, and Segal (2023)
% Load the CRSP and Compustat data needed to construct the variables;
load RECTRQ
load SALEQ
load SIC
load dates
load me

% Get rid of financials
RECTRQ(floor(SIC/1000)==6) = NaN;
RECTRQ(floor(SIC/10)==49) = NaN;

SALEQ(floor(SIC/1000)==6) = NaN;
SALEQ(floor(SIC/10)==49) = NaN;

me(floor(SIC/1000)==6) = NaN;
me(floor(SIC/10)==49) = NaN; 

% Get rid of zero sale firms;
SALET = SALEQ;
SALET(SALET==0) = NaN;
RS = (RECTRQ./SALET);
RS(RS==0) = NaN;
me(isnan(RS)) = NaN;
TC_Level = sum(RS.*me./sum(me, 2, 'omitnan'), 2, 'omitnan');
factors(5).factor = log(TC_Level)-log(lag(TC_Level, 1, nan));
factors(5).label = {'$\overline{R/S}$'};
clearvars -except factors dates

% IMC factor
load SIC
load ret
load dates
load me

% Download & read the classification file from Moto Yogo's website
% getGoogleDriveData('SIC_Durables.xls', '0BzR-ojpYuaFMQjFBaXo2RnlUTkU', 'Data\');
opts = detectImportOptions('SIC_Durables.xls', 'Sheet', 'SIC_Durables');
dataYogo = readtable('SIC_Durables.xls', opts);

% Classify the G and NX industries into consumption and investment:
indCons = ismember(dataYogo.Durability, {'G','NX'}) & ...
          ismember(dataYogo.SICCode, [780 7360 7370 7381 7382 8700 8730 8900 1240 2865]);
dataYogo.Durability(indCons) = repmat({'d.'}, sum(indCons), 1);
indInv = ismember(dataYogo.Durability, {'G','NX'}) & ...
         ~ismember(dataYogo.SICCode, [780 7360 7370 7381 7382 8700 8730 8900 1240 2865]);
dataYogo.Durability(indInv) = repmat({'I'}, sum(indInv), 1);

consumptionInvestment = zeros(size(SIC));

for i = 1:height(dataYogo)
    if ismember(dataYogo.Durability(i), {'d.','n.d.','s.'})
        consumptionInvestment(SIC==dataYogo.SICCode(i)) = 1;
    elseif ismember(dataYogo.Durability(i), {'I'})
        consumptionInvestment(SIC==dataYogo.SICCode(i)) = 2;
    end
end

res = runUnivSort(ret, consumptionInvestment, dates, me, 'plotFigure', 0, 'printResults', 0);
IMC = res.pret(:,end);
save Data\IMC IMC

factors(6).factor = IMC;
factors(6).label = {'IMC'};
clearvars -except factors dates

save Data\factors factors

%% Make industry-adjusted LPC

clear
clc

load FF17
load ret
load lpc

iaLPC = nan(size(ret));

for i = 1:17
    % Leave only the current industry
    tempLPC = lpc;
    tempLPC(FF17~=i) = nan;
    
    % Calculate the mean LPC for this industry in each month
    meanLPC = mean(tempLPC, 2, 'omitnan');
    rptdMeanLPC = repmat(meanLPC, 1, size(ret, 2));

    % Subtract the mean LPC 
    iaLPC(FF17==i) = lpc(FF17==i) - rptdMeanLPC(FF17==i);
end

% Store the variable
save Data\iaLPC iaLPC


%% Placebo LPC - [0, 50)

clear
clc

load ret
load nyse
load dates
load permno
load me
load paydex
load ff

% We'll start in 2005
s = find(dates==200501);

% Store a few constants
[nMonths, nStocks] = size(ret);

% initialize the LPC
lpc_0_50 = nan(nMonths, nStocks);

nroll = 60;

% Loop through the months
for i=s:nMonths 

    % Calculate LPC
    lpc_0_50(i,:) = sum(paydex(i-nroll+1:i, :) >= 0 & ...
                   paydex(i-nroll+1:i, :) < 50, 1) ./ ...
               sum(isfinite(paydex(i-nroll+1:i, :)));
end   

% Set observations for which we don't have paydex to NaN
lpc_0_50(isnan(paydex)) = nan;
    

save Data\lpc_0_50 lpc_0_50

%% Placebo LPC - [80, 100)

clear
clc

load ret
load nyse
load dates
load permno
load me
load paydex
load ff

% We'll start in 2005
s = find(dates==200501);

% Store a few constants
[nMonths, nStocks] = size(ret);

nroll = 60;

% initialize the LPC
lpc_80_100 = nan(nMonths, nStocks);

% Loop through the months
for i=s:nMonths 

    % Calculate LPC
    lpc_80_100(i,:) = sum(paydex(i-nroll+1:i, :) >= 80 & ...
                   paydex(i-nroll+1:i, :) < 100, 1) ./ ...
               sum(isfinite(paydex(i-nroll+1:i, :)));
end   

% Set observations for which we don't have paydex to NaN
lpc_80_100(isnan(paydex)) = nan;

save Data\lpc_80_100 lpc_80_100


%% Export data for Stata - Figure 2 and Tables 4 and 6

clear
clc

load z_paydex
load R
load IVOL
load me
load WhitedWu
load dates
load GP
load AT
load FinFirms
load permno
load lpc
load pmf
load market_share
load FQTR
load CAR3
load bm
load OIBDPQ
load ATQ
load AT
load FCF
load XRD

opCF = OIBDPQ./ATQ;
opCF(FinFirms==1)=nan;

fcf = FCF./AT;
fcf(FinFirms==1)=nan;

% Gross profitability
gp = GP./AT;
gp(FinFirms==1) = nan;


% Read in the ROE and distress signals from the novyMarxVelikovAnomalies.csv
[anoms, labels] = getAnomalySignals('novyMarxVelikovAnomalies.csv', 'permno', 'dates', ...
                                                'anomalyNames', {'roe','distress'});


% ROE 
r = strcmp(labels, 'roe');
roe = anoms(:, :, r);
roe(FinFirms==1) = nan;

% Distress
r = strcmp(labels, 'distress');
distress = -anoms(:, :, r);

% Whited-Wu
WhitedWu(FinFirms==1) = nan;

% Scale LPC
lpc = 100*lpc;

% bm(bm<=0) = nan;
load BEQ
bm = BEQ./me;
bm(bm<=0) = nan;

% CAR3
CAR3(isnan(FQTR)) = nan;

rdi = FillMonths(XRD)./me;

% Get GSCPI
opts = detectImportOptions('gscpi.xlsx');
data = readtable('gscpi.xlsx', opts);
data.Var1 = datetime(data.Var1);
data.Var1 = 100*year(data.Var1) + month(data.Var1);
[~, ia, ib] = intersect(dates, data.Var1);
gscpi = nan(size(dates));
gscpi(ia) = data.Var2(ib);

% Create table & export
[nMonths, nStocks] = size(AT);
nObs = nMonths * nStocks;
data = [reshape(z_paydex, nObs, 1) ...
        reshape(lpc, nObs, 1) ...
        reshape(WhitedWu, nObs, 1) ...
        reshape(FillMonths(pmf), nObs, 1) ...
        reshape(FillMonths(market_share), nObs, 1) ...
        reshape(distress, nObs, 1) ...
        reshape(R, nObs, 1) ...
        reshape(IVOL, nObs, 1) ...
        reshape(log(me), nObs, 1) ...
        reshape(FillMonths(gp), nObs, 1) ...
        reshape(roe, nObs, 1) ...
        reshape(CAR3, nObs, 1) ...        
        reshape(repmat(gscpi, 1, nStocks), nObs, 1) ...
        reshape(log(bm), nObs, 1) ...        
        reshape(bm, nObs, 1) ...        
        reshape(rdi, nObs, 1) ...        
        reshape(FillMonths(opCF), nObs, 1) ...        
        reshape(FillMonths(fcf), nObs, 1) ...        
        reshape(repmat(dates, 1, nStocks), nObs, 1) ...
        reshape(repmat(permno', nMonths, 1), nObs, 1)];
data(isnan(data(:,1)) & isnan(data(:,2)), :) = [];
data = array2table(data, 'VariableNames', {'zpaydex','lpc','WhitedWu','pmf','marketshare','distress','R','IVOL','logME','GP','roe','CAR3','GSCPI','logBM','BM','rdi','opCF','FCF','dates','permno'});

writetable(data,'Data\stata_data.csv');
