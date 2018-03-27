close all;
clear;

saveAggregatePlots = 0;
plotIndividualResults = 1;

%penetration = 0.2:0.05:0.7;
%plot_run = 'run29';
%plot_cases = {'run26_E1_PQ_P065','run26_E1_PQ_P07'};
plot_cases = {'run32_E1_PQ_P095'};

if exist('plot_run','var') % plot whole run, overries plot_cases
    
    D = dir(['data/' plot_run '*.mat']);
    N = length(D);

    for i=1:N
        
        caseNames{i} =  strrep( strrep( strrep( D(i).name,[plot_run '_'],'' ), '_', ' '), '.mat','' );
        
        % load data
        load(['data/' D(i).name]);
        
        f_values(i) = restab.f;
        outmsg{i} = Output.message;
        
        % Total curtailment (for each scenario and case in MW)
        curtail = table2array(restab.Curtail(:,setdiff( restab.Curtail.Properties.VariableNames,{'GEN','BUS'})));
        curtailedTotal = sum(curtail,1);
        curtailment_nom(i,:) = curtailedTotal;
        
        % Expected curtailment (for each case in MW)
        curtailment_nom_exp(i) = sum(table2array(restab.ExpCurtail(3,:))); % same as f-values
        
        % Curtailment (for each scenario and case in MW)
        wind = table2array(restab.Wind(:,setdiff(restab.Wind.Properties.VariableNames,{'GEN','BUS'})));
        windTotal = sum(wind,1);
        
        curtailment_pc(i,:) = 100*curtailedTotal./windTotal;
        
        % Expected curtailment (for each case in %)
        curtailment_pc_exp(i) = sum(curtailment_pc(i,:) .* table2array(restab.ExpCurtail(1,:)));
        
        % Transfer through corridors
        transferNorthSouth(i,:) = table2array(restab.transferFrom(1,setdiff( restab.transferFrom.Properties.VariableNames,{'CORRIDOR'})));
        
        % successful optimization
        success(i) = restab.success;     
    end
       
else % plot cases in plot_cases
    
    N = length(plot_cases);
    %caseNames = cell(N,1);
    
    for i=1:N
        
        D = dir(['data/' plot_cases{i} '.mat']);
        if ~isempty(D)
            
            load(['data/' D(1).name]);
            
            if plotIndividualResults
                plot_results(restab,optns);
            end
            
            caseNames{i} = strrep( strrep( D(1).name,[plot_cases{i} '_'],'' ), '_', ' ');
            
            f_values(i) = restab.f;
            outmsg{i} = Output.message;
            
            % Total curtailment (for each scenario and case in MW)
            curtail = table2array(restab.Curtail(:,setdiff( restab.Curtail.Properties.VariableNames,{'GEN','BUS'})));
            curtailedTotal = sum(curtail,1);
            curtailment_nom(i,:) = curtailedTotal;
            
            % Expected curtailment (for each case in MW)
            curtailment_nom_exp(i) = sum(table2array(restab.ExpCurtail(3,:))); % same as f-values
            
            % Curtailment (for each scenario and case in MW)
            wind = table2array(restab.Wind(:,setdiff(restab.Wind.Properties.VariableNames,{'GEN','BUS'})));
            windTotal = sum(wind,1);
            
            curtailment_pc(i,:) = 100*curtailedTotal./windTotal;
            
            % Expected curtailment (for each case in %)
            curtailment_pc_exp(i) = sum(curtailment_pc(i,:) .* table2array(restab.ExpCurtail(1,:)));
            
            % Transfer through corridors
            transferNorthSouth(i,:) = table2array(restab.transferFrom(1,setdiff( restab.transferFrom.Properties.VariableNames,{'CORRIDOR'})));
            
            % successful optimization
            success(i) = restab.success;
            
        else
            disp('No such file');
        end
    end
    
    plot_run = plot_cases{1};
    if contains( plot_run,'_' )
        plot_run = plot_run( 1:strfind( plot_run,'_' )-1 );
    end
end


%% PLOTS

% objective function
F1 = figure;
bar(categorical(caseNames),curtailment_nom_exp.');
ylabel('Curtailment (MW)');
xlabel('Case');
title(['Total expected curtailment' '\newline' strrep(plot_run,'_',' ')])


F2 = figure;
bar(1:size(curtailment_nom,2),curtailment_nom.');
xlabel('Contingency scenario');
L2 = legend(caseNames);
L2.Position = [0.25 0.25 0 0];
ylabel('Curtailment (MW)');
title(['Total curtailment' '\newline' strrep(plot_run,'_',' ')])

F3 = figure;
bar(1:size(curtailment_pc,2),curtailment_pc.');
ylabel('Curtailment (%)');
L3 = legend(caseNames);
L3.Position = [0.25 0.25 0 0];
xlabel('Contingency scenario');
title(['Total curtailment' '\newline' strrep(plot_run,'_',' ')]);

F4 = figure;
bar(categorical(caseNames),curtailment_pc_exp.');
ylabel('Curtailment (%)');
xlabel('Case');
title(['Total expected curtailment' '\newline' strrep(plot_run,'_',' ')])


F5 = figure;
bar(1:size(transferNorthSouth,2),transferNorthSouth.');
ylabel('Transfer (MW)');
L5 = legend(caseNames);
L5.Position = [0.25 0.25 0 0];
xlabel('Contingency scenario');
title(['Transfer from North to Central' '\newline' strrep(plot_run,'_',' ')]);
figs = [F1 F2 F3 F4 F5];

if exist('penetration','var')
    F6 = figure;
    plot(penetration,curtailment_nom_exp,'*');
    grid on;
    xlabel('Penetration (%)');
    ylabel('Curtailment (MW)');
    title(['Total expected curtailment' '\newline' strrep(plot_run,'_',' ')])
    figs = [figs F6];
    
end

F7 = figure;
bar(categorical(caseNames),success.');
xlabel('Case')
title(['Successful optimization' '\newline' strrep(plot_run,'_',' ')])

figs = [figs F7];
if exist('penetration','var')
    figNames = {'CurTotMW','CurMW','CurPC','CurTotPC','Transfer','CurTotMW_xy','Success'};
else
    figNames = {'CurTotMW','CurMW','CurPC','CurTotPC','Transfer','Success'};
end

if saveAggregatePlots
    for i=1:length(figs)
        thisFigureName = [plot_run '_' figNames{i}];
        saveas(figs(i),['figures/runs/' thisFigureName '.fig']);
        saveas(figs(i),['figures/runs/' thisFigureName '.png']);
    end
end
