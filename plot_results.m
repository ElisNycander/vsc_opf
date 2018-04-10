function [figureArray] = plot_results(table,optns)
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
n_legend_columns = 6;
%% make legends
ng = size(table.Pg,1);
genLabels = cell(1,ng);
for i=1:ng
   strIdx = num2str(table2array(table.Pg(i,'GEN')));
   %genLabels{i} = ['GEN ' strIdx];
   busStr = num2str(table2array(table.Pg(i,'BUS')));
   genLabels{i} =  ['GEN ' strIdx ' BUS ' busStr];
end

ncur = size(table.Beta,1);
curLabels = cell(1,ncur);
for i=1:ncur
    strIdx = num2str(table2array(table.Beta(i,'GEN')));
    busIdx = num2str(table2array(table.Beta(i,'BUS')));
    curLabels{i} = ['GEN ' strIdx ' BUS ' busIdx];
end

%% Pg for all generators
F1 = figure;
PgMat = table2array(table.Pg(:, setdiff(table.Pg.Properties.VariableNames,{'GEN','BUS','MIN','MAX'})));
PH1 = bar(1:N,PgMat');
% change distinguish wind farms from other generators
%for i=1:length(optns.gen.curtailableP)
%   PH1(optns.gen.curtailableP(i)).EdgeColor = [1 0 0];
%end
if ~isfield(optns,'plotGreyscale') || ~optns.plotGreyscale
    for i=1:length(PH1)
        
        % Edge: Black - variable, Red - curtailable, None - Fixed
        % determine type of generator
        if ismember( i,optns.gen.curtailableP ) % curtailable
            PH1(i).EdgeColor = [0.75 0 0];
        elseif ismember( i,optns.gen.fixedP ) % fixed
            PH1(i).EdgeColor = [0 0 0];
        else % variable
            PH1(i).EdgeColor = [1 1 1];
        end
        % determine location of generator
        % North: Blue, Central: Green, South: Yellow, External: Grey
        genbus = table2array( table.Qg( i,'BUS' ) );
        genbus2d = rem( genbus,100 );
        if genbus2d < 40 % North
            PH1(i).FaceColor = [0.2810 0.3228 0.9579];
        elseif genbus2d < 60 % Central
            PH1(i).FaceColor = [0.3406 0.8008 0.4789];
        elseif genbus2d < 70 % South
            PH1(i).FaceColor = [0.9610 0.8890 0.1537];
        else % External
            PH1(i).FaceColor = [0.8571 0.8571 0.8571];
        end
    end
    
else % greyscale color code
    
    generator_greyscale(PH1,optns,table);
end

grid on;
xlabel('Scenario');
ylabel('P (MW)')
columnlegend(n_legend_columns,genLabels,'Location','NorthWest');
%title(['Active power' '\newline' strrep(optns.caseName,'_',' ')])

ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

MaximizeFigureWindow();


%% Qg for all generators
F2 = figure;
PH2 = bar(1:N,table2array(table.Qg(:, ...
    setdiff(table.Qg.Properties.VariableNames,{'GEN','BUS','MIN','MAX'})))' ...
);
% for i=1:length(optns.gen.curtailableP)
%    PH2(optns.gen.curtailableP(i)).EdgeColor = [1 0 0];
% end
if ~isfield(optns,'plotGreyscale') || ~optns.plotGreyscale
    for i=1:length(PH2)
        % Edge: Black - variable, Red - curtailable, None - Fixed
        % determine type of generator
        if ismember( i,optns.gen.curtailableP ) % curtailable
            PH2(i).EdgeColor = [0.75 0 0];
        elseif ismember( i,optns.gen.fixedP ) % fixed
            PH2(i).EdgeColor = [0 0 0];
        else % variable
            PH2(i).EdgeColor = [1 1 1];
        end
        % determine location of generator
        % North: Blue, Central: Green, South: Yellow, External: Grey
        genbus = table2array( table.Qg( i,'BUS' ) );
        genbus2d = rem( genbus,100 );
        if genbus2d < 40 % North
            PH2(i).FaceColor = [0.2810 0.3228 0.9579];
        elseif genbus2d < 60 % Central
            PH2(i).FaceColor = [0.3406 0.8008 0.4789];
        elseif genbus2d < 70 % South
            PH2(i).FaceColor = [0.9610 0.8890 0.1537];
        else % External
            PH2(i).FaceColor = [0.8571 0.8571 0.8571];
        end
    end
else
    
    generator_greyscale(PH2,optns,table);
end

grid on;
xlabel('Scenario');
ylabel('Q (MW)')
columnlegend(n_legend_columns,genLabels,'Location','NorthWest');
%title(['Reactive power' '\newline' strrep(optns.caseName,'_',' ')])

MaximizeFigureWindow();

%% Wind Penetration
internalGens = false(size(table.Pg,1));
for i=1:length(internalGens)
    if ~ismember(table2array(table.Pg(i,'BUS')),optns.externalBuses)
        internalGens(i) = true;
    end
end


PgMat(isnan(PgMat)) = 0;
windPG = sum(PgMat(optns.gen.curtailableP,:),1);
totPG = sum(PgMat(internalGens,:),1);
penetration = windPG./totPG;

F3 = figure;
PH3 = bar(1:N,penetration*100);
grid on;
xlabel('Scenario');
ylabel('Wind Penetration (%)');
title(['Penetration' '\newline' strrep(optns.caseName,'_',' ')]);



%% Curtailment for Pcur
F4 = figure;
PH4 = bar(1:N, table2array(table.Curtail(:, ...
    setdiff(table.Curtail.Properties.VariableNames,{'GEN','BUS'})))' ...
);
grid on;
xlabel('Scenario');
ylabel('Curtailed energy (MW)');
legend(curLabels);
title(['Curtailment' '\newline' strrep(optns.caseName,'_',' ')])

%% Beta for Pcur
F5 = figure;
PH5 = bar(1:N, table2array(table.Beta(:, ...
    setdiff(table.Beta.Properties.VariableNames,{'GEN','BUS'})))' ...
);
grid on;
xlabel('Scenario');
ylabel('Beta');
legend(curLabels);
title(['Curtailment fraction beta' '\newline' strrep(optns.caseName,'_',' ')])

%% Total curtailment

curtail = table2array(table.Curtail(:,setdiff( table.Curtail.Properties.VariableNames,{'GEN','BUS'})));
curtailedTotal = sum(curtail,1);

F6 = figure;
PH6 = bar(1:N, curtailedTotal);
grid on;
xlabel('Scenario');
ylabel('Curtailment (MW)');
title(['Total curtailment' '\newline' strrep(optns.caseName,'_',' ')])

%% Total curtailment in %
wind = table2array(table.Wind(:,setdiff(table.Wind.Properties.VariableNames,{'GEN','BUS'})));
windTotal = sum(wind,1);

F7 = figure;
PH7 = bar(1:N, 100*curtailedTotal./windTotal);
grid on;
xlabel('Scenario');
ylabel('Curtailment (%)');
title(['Total curtailment' '\newline' strrep(optns.caseName,'_',' ')])


%% Plot transfers in corridors

transfers = table2array( table.transferFrom(:,setdiff( table.transferFrom.Properties.VariableNames,{'CORRIDOR'} )) );
F8 = figure;
PH8 = bar(1:N, transfers.');
transLables = {};
for i=1:size(transfers,1)
    transLables{i} = ['Corridor ' num2str(i)];
end
legend(transLables);
xlabel('Scenario');
ylabel('Active power transfer (MW)');
title(['Transfer through corridors' '\newline' strrep(optns.caseName,'_',' ')]);


%% Voltages

F9 = figure;
VmMat = table2array(table.Vm(:, setdiff(table.Vm.Properties.VariableNames,{'BUS','MIN','MAX'})));
PH9 = bar(1:N,VmMat');

% categories: 
% No generation, Wind, Conventional, Wind+Conventional: 1,2,3,4
catColors = [
    0 0 0; % black
    0 1 0; % green
    1 0 0; % red
    1 1 0 % yellow
];
vmLabels = {};
for i=1:size(VmMat,1)
    vmLabels{i} = ['Bus ' num2str( table2array( table.Vm(i,'BUS') ))];
end
legend(vmLabels);
columnlegend(n_legend_columns,vmLabels,'Location','NorthWest');

if ~isfield(optns,'plotGreyscale') || ~optns.plotGreyscale
    % set borders to mark categories
    for i=1:size(VmMat,1)
        bus = table2array( table.Vm(i,'BUS') );
        % find generators at this bus
        type = 1;
        for ii=1:size(table.Pg,1)
            if table2array(table.Pg(ii,'BUS')) == bus
                % check the type of this generator
                if ismember(ii,optns.gen.curtailableP) % wind
                    if type == 1
                        type = 2;
                    elseif type == 3
                        type = 4;
                    end
                else % conventional
                    if type == 1
                        type = 3;
                    elseif type == 2
                        type = 4;
                    end
                end
            end
        end
        %PH9(i).EdgeColor = catColors(type,:);
        PH9(i).FaceColor = catColors(type,:);
    end
    
else % greyscale
    voltage_greyscale(PH9,optns,table);
    
end

grid on;
xlabel('Scenario');
ylabel('V (pu)')

%title(['Voltage' '\newline' strrep(optns.caseName,'_',' ')])
ylim([0.85 1.15]);

MaximizeFigureWindow();


%% save figures
figureArray = [F1 F2 F3 F4 F5 F6 F7 F8 F9];
figNames = {'Pg','Qg','WindPenetration','CutailMW','CurtailFraction','CurtailTot','CurtailTotPc','Transfers','Vm'};
if optns.saveFigures
   for i=1:length(figureArray)
       thisFigureName = [optns.caseName '_' figNames{i}];
        saveas(figureArray(i),['figures/' thisFigureName '.fig']);
        saveas(figureArray(i),['figures/' thisFigureName '.png']);
   end
end

