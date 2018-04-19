%close all;
clear;

saveAggregatePlots = 0;
plotIndividualResults = 1;

runs = {'run41','run43','run42'};

run_legends = {'C1, base case', 'C2, no wind at 1021', 'C3, wind PF \geq 0.9'};
short_run_legends = {'C1','C2','C3'};
run_markers = {'.','*','o'};
run_marker_size = [6 4 4];
run_colors = [0 0 0;0 0 0;0 0 0];

penetration = 0.8:0.005:0.995;

curtail_pc_all = [];
windnet_all = [];
transfer_all = [];
curtail_all = [];
losses_all = [];


for j=1:length(runs)   

    plot_run = runs{j};
    
    D = dir(['data/' plot_run '*.mat']);
    N = length(D);
    
    caseNames = {};
    f_values = [];
    outmsg = {};
    curtailment_nom = [];
    curtailment_nom_exp = [];
    curtailment_pc = [];
    transferNorthSouth = [];
    losses = [];
    losses_pc = [];
    windnet = [];
    windscale = [];
    success = [];
    

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
        transferNorthSouth(i,:) = abs(table2array(restab.transferTo(1,setdiff( restab.transferTo.Properties.VariableNames,{'CORRIDOR'}))));
        
        
        cols = restab.Ploss.Properties.VariableNames;
        columns = {};
        for ii=1:length(cols)
            if contains(cols{ii},'MW')
                columns{end+1} = cols{ii};
            end
        end
        tmp = table2array( restab.Ploss(:,columns ), 1);
        tmp( isnan(tmp) ) = 0;
        losses(i,:) = sum( tmp, 1 );
        
        % sum total generation
        Pg = table2array( restab.Pg( :,setdiff( restab.Pg.Properties.VariableNames, {'GEN','BUS','MIN','MAX'} ) ) );
        Pg( isnan(Pg) ) = 0;
        ptot = sum( Pg, 1 );
        
        losses_pc(i,:) = losses(i,:) ./ ptot * 100;
        
        % net wind = wind - losses
        pwind = sum( Pg(optns.gen.curtailableP,:),1 );
        
        
        windnet(i,:) = pwind - losses(i,:);
        
        windscale(i) = pwind(1); % base case wind
        
        % successful optimization
        success(i) = restab.success;    
       
      
       

    end
    % save desired variables
    curtail_all = [curtail_all curtailment_nom(:,2)];
    curtail_pc_all = [curtail_pc_all curtailment_pc(:,2)];
    windnet_all = [windnet_all windnet(:,2)];
    transfer_all = [transfer_all transferNorthSouth];
    losses_all = [losses_all losses];
    
    
end

    % create legend
    run_legends_x_contingencies = {};
    ij = 1;
    for i=1:length(runs)
        N = size(restab.ExpCurtail,2);
        for j=1:N
            run_legends_x_contingencies{ij} = [short_run_legends{i} ', S' num2str(j)];
            ij = ij+1;
        end
    end

%% PLOTS,penetration,windscale)

F13 = figure;
PH13 = plot(penetration, curtail_all); % only high-wind scenarios
ylabel('Curtailment (MW)');
L13 = legend(run_legends);
grid on;
double_x_axis(F13,penetration,windscale);
change_markers(PH13,run_markers,run_marker_size);
%tight_axis(F13);

F14 = figure;
PH14 = plot(penetration,windnet_all);
grid on;
ylabel('P_{wind} - P_{loss} (MW)');
legend(run_legends);
grid on;
%axis tight
%ylim([3500 4500]);
double_x_axis(F14,penetration,windscale);
change_markers(PH14,run_markers,run_marker_size);
%tight_axis(F14);

F15 = figure;
PH15 = plot(penetration,transfer_all);
grid on;
ylabel('Transfer (MW)');
legend(run_legends_x_contingencies);
double_x_axis(F15,penetration,windscale);
change_markers(PH15,run_markers,run_marker_size);

F16 = figure;
PH16 = plot(penetration,losses_all);
grid on;
ylabel('Losses (MW)');
legend(run_legends_x_contingencies);
double_x_axis(F16,penetration,windscale);
change_markers(PH16,run_markers,run_marker_size);
% if saveAggregatePlots
%     for i=1:length(figs)
%         thisFigureName = [plot_run '_' figNames{i}];
%         saveas(figs(i),['figures/runs/' thisFigureName '.fig']);
%         saveas(figs(i),['figures/runs/' thisFigureName '.png']);
%     end
% end
