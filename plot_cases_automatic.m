close all;
clear;

saveAggregatePlots = 0;
plotIndividualResults = 1;
grey_plot = 0;

%penetration = 0.8:0.005:0.995;
%plot_run = 'run38';
%plot_cases = {'run26_E1_PQ_P065','run26_E1_PQ_P07'};
plot_cases = {'run43_E3_PQ_P0995'};

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
        transferNorthSouth(i,:) = table2array(restab.transferTo(1,setdiff( restab.transferFrom.Properties.VariableNames,{'CORRIDOR'})));
        
        
        cols = restab.Ploss.Properties.VariableNames;
        columns = {};
        for ii=1:length(cols)
            if contains(cols{ii},'MW')
                columns{end+1} = cols{ii};
            end
        end
        losses(i,:) = sum( table2array( restab.Ploss(:,columns) ), 1 );
        
        % sum total generation
        Pg = table2array( restab.Pg( :,setdiff( restab.Pg.Properties.VariableNames, {'GEN','BUS','MIN','MAX'} ) ) );
        ptot = sum( Pg, 1 );
        
        losses_pc(i,:) = losses(i,:) ./ ptot * 100;
        
        % net wind = wind - losses
        pwind = sum( Pg(optns.gen.curtailableP,:),1 );
        
        
        windnet(i,:) = pwind - losses(i,:);
        
        windscale(i) = pwind(1); % base case wind
        
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
            optns.plotGreyscale = grey_plot;
            
            if plotIndividualResults
                figArray = plot_results_article(restab,optns);
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
            transferNorthSouth(i,:) = table2array(restab.transferTo(1,setdiff( restab.transferFrom.Properties.VariableNames,{'CORRIDOR'})));
            
            
            cols = restab.Ploss.Properties.VariableNames;
            columns = {};
            for ii=1:length(cols)
                if contains(cols{ii},'MW')
                    columns{end+1} = cols{ii};
                end
            end
            losses(i,:) = sum( table2array( restab.Ploss(:,columns) ), 1 );
            
            % sum total generation
            Pg = table2array( restab.Pg( :,setdiff( restab.Pg.Properties.VariableNames, {'GEN','BUS','MIN','MAX'} ) ) );
            ptot = sum( Pg, 1 );
            
            losses_pc(i,:) = losses(i,:) ./ ptot * 100;
            
            % net wind = wind - losses
            pwind = sum( Pg(optns.gen.curtailableP,:),1 );
            
            windnet(i,:) = pwind - losses(i,:);
            
            windscale(i) = pwind(1); % base case wind
            
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

contingencyLabls = {};
for i=1:size(restab.ExpCurtail,2)
    contingencyLabls{i} = ['Contingency ' num2str(i)];
end

%% PLOTS
% 
% % objective function
% F1 = figure;
% bar(categorical(caseNames),curtailment_nom_exp.');
% ylabel('Curtailment (MW)');
% xlabel('Case');
% title(['Total expected curtailment' '\newline' strrep(plot_run,'_',' ')])
% 
% 
% F2 = figure;
% bar(1:size(curtailment_nom,2),curtailment_nom.');
% xlabel('Contingency scenario');
% L2 = legend(caseNames);
% L2.Position = [0.25 0.25 0 0];
% ylabel('Curtailment (MW)');
% title(['Total curtailment' '\newline' strrep(plot_run,'_',' ')])
% 
% F3 = figure;
% bar(1:size(curtailment_pc,2),curtailment_pc.');
% ylabel('Curtailment (%)');
% L3 = legend(caseNames);
% L3.Position = [0.25 0.25 0 0];
% xlabel('Contingency scenario');
% title(['Total curtailment' '\newline' strrep(plot_run,'_',' ')]);
% 
% F4 = figure;
% bar(categorical(caseNames),curtailment_pc_exp.');
% ylabel('Curtailment (%)');
% xlabel('Case');
% title(['Total expected curtailment' '\newline' strrep(plot_run,'_',' ')])
% 
% 
% F5 = figure;
% bar(1:size(transferNorthSouth,2),transferNorthSouth.');
% ylabel('Transfer (MW)');
% L5 = legend(caseNames);
% L5.Position = [0.25 0.25 0 0];
% xlabel('Contingency scenario');
% title(['Transfer from North to Central' '\newline' strrep(plot_run,'_',' ')]);
% figs = [F1 F2 F3 F4 F5];
% 
% if exist('penetration','var')
%     F6 = figure;
%     plot(penetration,curtailment_nom_exp,'*');
%     grid on;
%     xlabel('Penetration (%)');
%     ylabel('Curtailment (MW)');
%     title(['Total expected curtailment' '\newline' strrep(plot_run,'_',' ')])
%     figs = [figs F6];
%     
% end
% 
% F7 = figure;
% bar(categorical(caseNames),success.');
% xlabel('Case')
% title(['Successful optimization' '\newline' strrep(plot_run,'_',' ')])
% 
% figs = [figs F7];
% if exist('penetration','var')
%     figNames = {'CurTotMW','CurMW','CurPC','CurTotPC','Transfer','CurTotMW_xy','Success'};
% else
%     figNames = {'CurTotMW','CurMW','CurPC','CurTotPC','Transfer','Success'};
% end
% 
% % bar(1:size(curtailment_nom,2),curtailment_nom.');
% % xlabel('Contingency scenario');
% % L2 = legend(caseNames);
% % L2.Position = [0.25 0.25 0 0];
% % ylabel('Curtailment (MW)');
% % title(['Total curtailment' '\newline' strrep(plot_run,'_',' ')])
% 
% F8 = figure;
% bar(1:size(losses,2),losses.');
% xlabel('Contingency scenario');
% ylabel('MW');
% title(['Active power losses' '\newline' strrep(plot_run,'_',' ')]);
% 
% F9 = figure;
% bar(1:size(losses_pc,2),losses_pc.');
% xlabel('Contingency scenario');
% ylabel('%');
% title(['Active power losses' '\newline' strrep(plot_run,'_',' ')]);
% 
% figs = [figs F8 F9];
% figNames{end+1} = 'LossMW';
% figNames{end+1} = 'LossPC';
% 
% if exist('penetration','var') % plot losses as line plots
%    
%     
%     F10 = figure;
%     plot(penetration,losses);
%     grid on;
%     xlabel('Penetration (%)')
%     ylabel('Losses (MW)')
%     legend(contingencyLabls);
%     
%     F11 = figure; 
%     plot(penetration,losses_pc);
%     grid on
%     xlabel('Penetration (%)');
%     ylabel('Losses (MW)');
%     legend(contingencyLabls);
% 
%     fig = [figs F10 F11];
%     figNames{end+1} = 'LossXYMW';
%     figNames{end+1} = 'LossXYPC';
%     
%     F12 = figure;
%     plot(penetration,windnet);
%     grid on;
%     ylabel('P_{wind} - P_{loss} (MW)');
%     legend(contingencyLabls);
%     title('Net wind contribution');
%     double_x_axis(F12,penetration,windscale)
% end
% 
% 
% 
% if saveAggregatePlots
%     for i=1:length(figs)
%         thisFigureName = [plot_run '_' figNames{i}];
%         saveas(figs(i),['figures/runs/' thisFigureName '.fig']);
%         saveas(figs(i),['figures/runs/' thisFigureName '.png']);
%     end
% end

% change font size
% for i=1:length(figArray)
%    
%    figArray(i).CurrentAxes.XLabel.FontSize = 14;
%    figArray(i).CurrentAxes.YLabel.FontSize = 14;
%    figArray(i);
%    xt = get(gca,'XTick');
%    set(gca,'FontSize',14);
%    yt = get(gca,'YTick');
%    set(gca,'FontSize',14);
%    % figs(i).CurrentAxes.XLabel.FontSize = 14; 
% end
