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
N = size(table.Vm,2)-1;

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
figure;
PH1 = bar(1:N,table2array(table.Pg(:, ...
    setdiff(table.Pg.Properties.VariableNames,{'GEN','BUS'})))' ...
);
grid on;
xlabel('Contingency scenario');
ylabel('Pg (MW)')
legend(genLabels);
title('Active power')

%% Qg for all generators
figure;
PH2 = bar(1:N,table2array(table.Qg(:, ...
    setdiff(table.Qg.Properties.VariableNames,{'GEN','BUS'})))' ...
);
grid on;
xlabel('Contingency scenario');
ylabel('Qg (MW)')
legend(genLabels);
title('Reactive power')

%% Curtailment for Pcur
figure;
PH3 = bar(1:N, table2array(table.Curtail(:, ...
    setdiff(table.Curtail.Properties.VariableNames,{'GEN','BUS'})))' ...
);
grid on;
xlabel('Contingency scenario');
ylabel('Curtailed energy (MW)');
legend(curLabels);
title('Curtailment')

%% Beta for Pcur
figure;
PH4 = bar(1:N, table2array(table.Beta(:, ...
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

figure;
PH5 = bar(1:N, curtailedTotal);
grid on;
xlabel('Contingency scenario');
ylabel('Curtailment (MW)');
title('Total curtailment')

%% Total curtailment in %
wind = table2array(table.Wind(:,setdiff(table.Wind.Properties.VariableNames,{'GEN','BUS'})));
windTotal = sum(wind,1);

figure;
PH6 = bar(1:N, 100*curtailedTotal./windTotal);
grid on;
xlabel('Contingency scenario');
ylabel('Curtailment (%)');
title('Total curtailment')

