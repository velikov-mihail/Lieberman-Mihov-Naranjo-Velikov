function printFullSamplePerformanceTable(var, beg, numptf, nyseInd)

load dates
load ret
load me
load nyse
load ff

if nyseInd~=0
    ind = makeUnivSortInd(var, numptf, NYSE);
else
    ind = makeUnivSortInd(var, numptf, NYSE);
end

s = find(dates==beg);

res(1).res = runUnivSort(ret, ind, dates, me, 'timePeriod', beg, ...
                                              'plotFigure', 0, ...
                                              'printResults', 0, ...
                                              'factorModel', 1);
res(2).res = runUnivSort(ret, ind, dates, me, 'timePeriod', beg, ...
                                              'plotFigure', 0, ...
                                              'printResults', 0, ...
                                              'factorModel', 3);
res(3).res = runUnivSort(ret, ind, dates, me, 'timePeriod', beg, ...
                                              'plotFigure', 0, ...
                                              'printResults', 0, ...
                                              'factorModel', 4);
res(4).res = runUnivSort(ret, ind, dates, me, 'timePeriod', beg, ...
                                              'plotFigure', 0, ...
                                              'printResults', 0, ...
                                              'factorModel', 5);
res(5).res = runUnivSort(ret, ind, dates, me, 'timePeriod', beg, ...
                                              'plotFigure', 0, ...
                                              'printResults', 0, ...
                                              'factorModel', 6);

% clc

heads=[{'$r^e$'},{'$\alpha^{\text{CAPM}}$'},{'$\alpha^{\text{FF3}}$'},{'$\alpha^{\text{FF3+UMD}}$'},{'$\alpha^{\text{FF5}}$'},{'$\alpha^{\text{FF6}}$'}];


a = [res(1).res.xret res(1).res.alpha res(2).res.alpha res(3).res.alpha res(4).res.alpha res(5).res.alpha]';
tA = [res(1).res.txret res(1).res.talpha res(2).res.talpha res(3).res.talpha res(4).res.talpha res(5).res.talpha]';


% a(:,end)=-a(:,end);
% tA(:,end)=-tA(:,end);

mat2Tex(a, tA, heads, 2);


fprintf('\n\n');


heads = [{'$\beta_{\text{MKT}}$'},{'$\beta_{\text{SMB}}$'},{'$\beta_{\text{HML}}$'},{'$\beta_{\text{RMW}}$'},{'$\beta_{\text{CMA}}$'},{'$\beta_{\text{UMD}}$'}];
a = [res(5).res.factorLoadings(1).b res(5).res.factorLoadings(2).b res(5).res.factorLoadings(3).b res(5).res.factorLoadings(4).b res(5).res.factorLoadings(5).b res(5).res.factorLoadings(6).b]';
tA = [res(5).res.factorLoadings(1).t res(5).res.factorLoadings(2).t res(5).res.factorLoadings(3).t res(5).res.factorLoadings(4).t res(5).res.factorLoadings(5).t res(5).res.factorLoadings(6).t]';

% a(:,end)=-a(:,end);
% tA(:,end)=-tA(:,end);

mat2Tex(a,tA,heads,2);
