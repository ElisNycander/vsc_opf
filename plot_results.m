function [PH1,PH2,PH3,PH4,PH5,PH6] = plot_results(table,optns)
% close all;
% clear;
% load('plotData.mat');

%% options
%plotGens = [1:6];

%% unpack
%mpc = get_mpc(om);
%nCont = mpc.contingencies.nContingencies;
%nScen = mpc.contingencies.nScenarios;
%N = mpc.contingencies.N;
N = size(table.Vm,2)-3;

%% make legends
ng = size(table.Pg,1);
genLabels = cell(1,ng);
for i=1:ng
   strIdx = num2str(table2array(table.Pg(i,'GEN')));
   genLabels{i} = ['GEN ' strIdx];
end

ncur = size(table.Beta,1);
curLabels = cell(1,ncur);
for i=1:ncur
    strIdx = num2str(table2array(table.Beta(i,'GEN')));
    curLabels{i} = ['GEN ' strIdx];
end

%% Pg for all generators
F1 = figure;
PgMat = table2array(table.Pg(:, setdiff(table.Pg.Properties.VariableNames,{'GEN','BUS','MIN','MAX'})));
PH1 = bar(1:N,PgMat');
% change distinguish wind farms from other generators
for i=1:length(optns.gen.curtailableP)
   PH1(optns.gen.curtailableP(i)).EdgeColor = [1 0 0]; 
end

grid on;
xlabel('Contingency scenario');
ylabel('Pg (MW)')
legend(genLabels);
title('Active power')

MaximizeFigureWindow();


%% Qg for all generators
F2 = figure;
PH2 = bar(1:N,table2array(table.Qg(:, ...
    setdiff(table.Qg.Properties.VariableNames,{'GEN','BUS','MIN','MAX'})))' ...
);
for i=1:length(optns.gen.curtailableP)
   PH2(optns.gen.curtailableP(i)).EdgeColor = [1 0 0]; 
end

grid on;
xlabel('Contingency scenario');
ylabel('Qg (MW)')
legend(genLabels);
title('Reactive power')

MaximizeFigureWindow();

%% Wind Penetration
PgMat(isnan(PgMat)) = 0;
windPG = sum(PgMat(optns.gen.curtailableP,:),1);
totPG = sum(PgMat(),1);
penetration = windPG./totPG;

F3 = figure;
PH3 = bar(1:N,penetration*100);
grid on;
xlabel('Scenario');
ylabel('Wind Penetration (%)');
title('Penetration');



%% Curtailment for Pcur
F4 = figure;
PH4 = bar(1:N, table2array(table.Curtail(:, ...
    setdiff(table.Curtail.Properties.VariableNames,{'GEN','BUS'})))' ...
);
grid on;
xlabel('Contingency scenario');
ylabel('Curtailed energy (MW)');
legend(curLabels);
title('Curtailment')

%% Beta for Pcur
F5 = figure;
PH5 = bar(1:N, table2array(table.Beta(:, ...
    setdiff(table.Beta.Properties.VariableNames,{'GEN','BUS'})))' ...
);
grid on;
xlabel('Contingency scenario');
ylabel('Beta');
legend(curLabels);
title('Curtailment fraction beta')

%% Total curtailment

curtail = table2array(table.Curtail(:,setdiff( table.Curtail.Properties.VariableNames,{'GEN','BUS'})));
curtailedTotal = sum(curtail,1);

F6 = figure;
PH6 = bar(1:N, curtailedTotal);
grid on;
xlabel('Contingency scenario');
ylabel('Curtailment (MW)');
title('Total curtailment')

%% Total curtailment in %
wind = table2array(table.Wind(:,setdiff(table.Wind.Properties.VariableNames,{'GEN','BUS'})));
windTotal = sum(wind,1);

F7 = figure;
PH7 = bar(1:N, 100*curtailedTotal./windTotal);
grid on;
xlabel('Contingency scenario');
ylabel('Curtailment (%)');
title('Total curtailment')

figureArray = [F1 F2 F3 F4 F5 F6 F7];
figNames = {'Pg','Qg','WindPenetration','CutailMW','CurtailFraction','CurtailTot','CurtailTotPc'};
if optns.saveFigures
   for i=1:length(figureArray)
       thisFigureName = [optns.caseName '_' figNames{i}];
        saveas(figureArray(i),['figures/' thisFigureName '.fig']);
        saveas(figureArray(i),['figures/' thisFigureName '.png']);
   end
end

