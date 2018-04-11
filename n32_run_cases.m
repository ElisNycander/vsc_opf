close all;
clear;

rerun_simulations = 1;
base_identifier = 'run44';

windScenario = 'E3';
penetration = 0.98:0.005:0.995;
%PQ_constraints = [0 1];
PQ_constraints = 1;
Q_wind = 0;
save_plots = 0;
show_plots = 0;

%% run simulations
jkl = 1;
for j=1:length(PQ_constraints)
    for k=1:length(Q_wind)
        
        for l=1:length(penetration)
            
            caseNames{jkl} = windScenario;
            if PQ_constraints(j) == 1
                caseNames{jkl} = [caseNames{jkl} '_PQ'];
            end
            if Q_wind(k) == 1
                caseNames{jkl} = [caseNames{jkl} '_Qwind'];
            end
            % add
            caseNames{jkl} = [caseNames{jkl} '_P' strrep( num2str(penetration(l) ),'.','')];
            caseNames{jkl} = [base_identifier '_' caseNames{jkl}];
            
            if rerun_simulations
                run_opf_n32(caseNames{jkl},windScenario,PQ_constraints(j),Q_wind(k),penetration(l),save_plots);
            end
            
            jkl = jkl + 1;
            
        end
    end
end

success = [];
f_values = [];
outmsg = {};
curtailment_nom = [];
curtailment_nom_exp = zeros(length(caseNames),1); % should be equal to f value
curtailment_pc = [];
curtailment_pc_exp = zeros(length(caseNames),1);
%% collect results
for i=1:length(caseNames)
    
    % load data
    load(['data/' caseNames{i} '.mat']);
    
    if show_plots
       plot_results(restab,optns); 
    end
    
    f_values(i) = restab.f;
    outmsg{i} = Output.message;
    
    %% Total curtailment (for each scenario and case in MW)
    curtail = table2array(restab.Curtail(:,setdiff( restab.Curtail.Properties.VariableNames,{'GEN','BUS'})));
    curtailedTotal = sum(curtail,1);
    curtailment_nom(i,:) = curtailedTotal;

    %% Expected curtailment (for each case in MW)
    curtailment_nom_exp(i) = sum(table2array(restab.ExpCurtail(3,:))); % same as f-values
    
    %% Curtailment (for each scenario and case in MW)
    wind = table2array(restab.Wind(:,setdiff(restab.Wind.Properties.VariableNames,{'GEN','BUS'})));
    windTotal = sum(wind,1);
    
    curtailment_pc(i,:) = 100*curtailedTotal./windTotal;
    
    %% Expected curtailment (for each case in %)
    curtailment_pc_exp(i) = sum(curtailment_pc(i,:) .* table2array(restab.ExpCurtail(1,:)));
    
    %% Transfer through corridors
    transferNorthSouth(i,:) = table2array(restab.transferFrom(1,setdiff( restab.transferFrom.Properties.VariableNames,{'CORRIDOR'})));
    
    %% successful optimization
    success(i) = restab.success;
end

caseNames = strrep(caseNames,'_',' ');

% objective function
F1 = figure;
bar(categorical(caseNames),curtailment_nom_exp.');
ylabel('Curtailment (MW)');
xlabel('Case');
title(['Total expected curtailment' '\newline' strrep(base_identifier,'_',' ')])



F2 = figure;
bar(1:size(curtailment_nom,2),curtailment_nom.');
xlabel('Contingency scenario');
L2 = legend(caseNames);
L2.Position = [0.25 0.25 0 0];
ylabel('Curtailment (MW)');
title(['Total curtailment' '\newline' strrep(base_identifier,'_',' ')])

F3 = figure;
bar(1:size(curtailment_pc,2),curtailment_pc.');
ylabel('Curtailment (%)');
L3 = legend(caseNames);
L3.Position = [0.25 0.25 0 0];
xlabel('Contingency scenario');
title(['Total curtailment' '\newline' strrep(base_identifier,'_',' ')]);

F4 = figure;
bar(categorical(caseNames),curtailment_pc_exp.');
ylabel('Curtailment (%)');
xlabel('Case');
title(['Total expected curtailment' '\newline' strrep(base_identifier,'_',' ')])


F5 = figure;
bar(1:size(transferNorthSouth,2),transferNorthSouth.');
ylabel('Transfer (MW)');
L5 = legend(caseNames);
L5.Position = [0.25 0.25 0 0];
xlabel('Contingency scenario');
title(['Transfer from North to Central' '\newline' strrep(base_identifier,'_',' ')]);
figs = [F1 F2 F3 F4 F5];

if exist('penetration','var')
    F6 = figure;
    plot(penetration,curtailment_nom_exp,'*');
    grid on;
    xlabel('Penetration (%)');
    ylabel('Curtailment (MW)');
    title(['Total expected curtailment' '\newline' strrep(base_identifier,'_',' ')])
    figs = [figs F6];
    
end

F7 = figure;
bar(categorical(caseNames),success.');
xlabel('Case')
title(['Successful optimization' '\newline' strrep(base_identifier,'_',' ')])

% figs = [figs F7];
% if exist('penetration','var')
%     figNames = {'CurTotMW','CurMW','CurPC','CurTotPC','Transfer','CurTotMW_xy','Success'};
% else
%     figNames = {'CurTotMW','CurMW','CurPC','CurTotPC','Transfer','Success'};
% end
% 
% for i=1:length(figs)
%     thisFigureName = [base_identifier '_' figNames{i}];
%     saveas(figs(i),['figures/runs/' thisFigureName '.fig']);
%     saveas(figs(i),['figures/runs/' thisFigureName '.png']);
% end

