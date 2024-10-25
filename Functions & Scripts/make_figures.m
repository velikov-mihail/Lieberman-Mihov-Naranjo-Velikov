%% Make Figures path

clear
clc

figuresPath = [pwd, filesep, 'Figures'];

if ~exist(figuresPath, 'dir')
  mkdir(figuresPath);
  addpath(genpath(figuresPath));
end


%% Figure 1: Liquid Holdings Group, Inc.

clear
clc

load ret
load nyse
load dates
load me
load paydex
load ff
load dret
load ddates
load RDQ
load dprc
load permno

% Store a few constants
axisSize = 25;
legendSize = 30;
lineWidth = 3;
markerSize = 10;
lightGray = [0.66,0.66,0.66];

% Create the datetime for the x axis
pddates = datetime(ddates, 'ConvertFrom', 'yyyyMMdd');

% Start the figure
h = figure('visible','off');

% Set the font
set(0,'DefaultAxesFontName', 'Times New Roman')

% Find the company and the dates
c = find(permno==14053);
r = find(ddates<20131201, 1, 'last');

% Find the 6-month range around end of 2013 range 
s = max( find(pddates>=(pddates(r)-calmonths(3)), 1, 'first'), ...
         find(isfinite(dprc(:,c)), 1, 'first'));
e = find(pddates<=(pddates(r)+calmonths(3)), 1, 'last');

% Create the two series we'll plot
x = ddates(s:e);
y1 = abs(dprc(s:e,c));
y2 = nan(size(y1));
for i=1:length(y2)
    p = find(dates<floor(x(i)/100), 1, 'last');
    y2(i) = paydex(p,c);
end
% Plot them 
x = pddates(s:e);

% Left axis
yyaxis left
plot(x, y1, 'LineStyle', '--', ....
            'color', lightGray, ...
            'LineWidth', lineWidth, ...
            'MarkerSize', markerSize);
ylabel('LIQD price ($)');

% Right axis
yyaxis right
plot(x, y2, 'LineStyle', '-', ...
            'color', 'k', ...
            'LineWidth', lineWidth, ...
            'MarkerSize', markerSize);
ylabel('LIQD paydex score');

% legend 
legend('LIQD Price ($)','LIQD Paydex Score','FontSize',legendSize);

% Clean up 
ylim([0 100]);
hold on;
SP = pddates(find(pddates<=datetime(2013,11,30), 1, 'last')) + caldays(2);
yy = get(gca, 'YLim');
xx = get(gca, 'XLim');
set(gca, 'XLim', [datetime(2013,9,1) datetime(2014,3,1)]);
r2 = pddates( find(pddates >= datetime(RDQ(find(RDQ(:,c) > ddates(r), 1, 'first'), c), 'ConvertFrom', 'yyyyMMdd'), 1, 'first'));

% Add the text
x = [0.27 0.345] + 0.155;
y = [0.4 0.66];
message = sprintf('LIQD enters short\n portfolio on 11/29/2013');
annotation('textarrow', x, y, 'String', message, 'FontSize', legendSize*.8)

% SEt the line width
hline = flipud(findobj(h, 'type', 'line'));
set(hline(1),'LineWidth', lineWidth);
set(hline(2),'LineWidth', lineWidth); 

% Final touch-ups
ax = gca;
ax.YAxis(1).Color = lightGray;
ax.YAxis(2).Color = 'k';
set(gca,'FontSize',axisSize');

set(gcf, 'PaperPositionMode', 'auto');
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
set(gca,'LooseInset',get(gca,'TightInset'))
set(gcf,'units','inches')
orient(gcf,'landscape')
set(gcf,'PaperSize',[20 10]); %set the paper size to what you want  

% Save the figure
print('Figures/liqdFigure','-dpdf','-fillpage')

%% Figure 3: Target & Walmart paydex

clear
clc

load paydex
load permno
load dates
load pdates

% Store a few constants
axisSize = 25;
lineWidth = 3;
markerSize = 10;
legendSize = 30;
darkGray = [0.33,0.33,0.33];

% Find the two companies' permnos
c1 = find(permno==49154);
c2 = find(permno==55976);

% Find the dates range
s = find(dates==200601);
e = find(dates==201912);

% Plot the figure
h = figure('visible','off');

% Set the font
set(h, 'DefaultAxesFontName', 'Times New Roman')

x = pdates;
plot(x(s:e), paydex(s:e,c1), 'LineStyle', ':', ...
                             'color', darkGray, ...
                             'LineWidth', lineWidth, ...
                             'MarkerSize', markerSize);
hold on;
plot(x(s:e), paydex(s:e,c2), 'LineStyle', '-', ...
                             'color', 'k', ...
                             'LineWidth', lineWidth, ...
                             'MarkerSize', markerSize);
ylabel('Paydex score');
legend('Target','Walmart','FontSize',legendSize);
ylim([70 80]);
hold on;

% Final touch-ups
hline = flipud(findobj(h, 'type', 'line'));
set(hline(1),'LineWidth', lineWidth);
set(hline(2),'LineWidth', lineWidth); 
set(gca,'FontSize',axisSize');
set(gcf, 'PaperPositionMode', 'auto');
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
set(gca,'LooseInset',get(gca,'TightInset'))
set(gcf,'units','inches')
orient(gcf,'landscape')
set(gcf,'PaperSize',[20 10]); %set the paper size to what you want  

% Store the figure
print('Figures/mktPowerFigure','-dpdf','-fillpage')

%% Figure 4: Paydex vs market cap

clear
clc

load me
load paydex

% Store a few constants
axisSize = 25;
titleSize=30;
markerSize = 50;
[nMonths, nStocks] = size(me);
nObs = nMonths * nStocks;

% Plot the figure
h = figure('visible','off');

% Set the font
set(h, 'DefaultAxesFontName', 'Times New Roman')

% Reshape & drop NaNs
rshpdData = array2table([reshape(me,nObs,1) reshape(paydex,nObs,1)], 'VariableNames', {'me','paydex'});
indToDrop = isnan(rshpdData.me) | isnan(rshpdData.paydex);
rshpdData(indToDrop, :) = [];

% Plot
scatter(rshpdData.paydex,rshpdData.me/1000, markerSize, 'filled', 'k');
xlabel('Paydex', 'FontSize', axisSize);
ylabel('Market capitalization ($ bln)', 'FontSize', axisSize);
title('Paydex vs market capitalization', 'FontSize', titleSize ,...
                                         'FontWeight', 'normal');

% Final touch-ups
set(gcf, 'PaperPositionMode', 'auto');
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
set(gca,'LooseInset',get(gca,'TightInset'))
set(gcf,'units','inches')
set(gca,'FontSize',axisSize');
orient(gcf,'landscape')
set(gcf,'PaperSize',[20 10]); %set the paper size to what you want  

% Store
print('Figures/mcapPaydexScatter','-dpdf','-fillpage')

