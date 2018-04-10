function [figureArray] = plot_results(table,optns)
% close all;
% clear;
% load('plotData.mat');

%% options
%plotGens = [1:6];

N = size(table.Vm,2)-3;
n_legend_columns = 6;

% default values
if ~isfield(optns.plot,'plots')
    plots = struct();
    [plots.Pg,plots.Qg,plots.Vm,plots.Curtail, ...
        plots.CurtailTotal,plots.CurtailTotalPercent,plots.Beta,...
        plots.Transfer,plots.Penetration] = deal(...
        1, 1, 1, 0, ...
        0, 0, 0, ...
        0, 0 );
else
    plots = optns.plot.plots;
end

if ~isfield(optns.plot,'plot_titles')
    plot_titles = 1;
else
    plot_titles = optns.plot.plot_titles;
end

if ~isfield(optns,'scenario_names')
    barx = categorical(1:N);
    %scenario_names = categorical(1:N)
else
   assert(N == length(optns.scenario_names),'''optns.scenario_names'' must be of same length as number of scenarios');
   scenario_names = categorical(optns.scenario_names);
   barx = reordercats(scenario_names,optns.scenario_names);
end

figureArray = [];


%% unpack
%mpc = get_mpc(om);
%nCont = mpc.contingencies.nContingencies;
%nScen = mpc.contingencies.nScenarios;
%N = mpc.contingencies.N;

yoffset = 0.01; % offsets for position of text above bars
xoffset = 0.0;

fwidth = 1250;
fheight = 305;
tx_margin = 1e-9;
height_adjust = 15; % increase window size to prevent text boxes from overlapping
% greyscale colors for different areas
%color_factor = linspace(0.2,0.9,4);
%color_factor = linspace(0.9,0.2,4);
color_factor = linspace(1,0.3,4);
%color_factor = ones(1,4);
busFontSize = 10;


ann_strs2 = {'Variable','Fixed'};
edge_styles = {'-','none'};
ann_strs = {'North','Central','Southwest','External'};
legend_greyfactor = 0.8;
% colors = [
%     0.2810 0.3228 0.9579
%     0.3406 0.8008 0.4789
%     0.9610 0.8890 0.1537
%     0.8571 0.8571 0.8571
% ];
%
% colors = [
%     0.3 0.3228 1
%     0.3406 0.8008 0.4789
%     0.9610 0.8890 0.1537
%     0.8571 0.8571 0.8571
% ];

% colors = [ % different strength
% 	213 255 255
% 	87  235 153
% 	198 176  42
% 	91   91  91
% ]/255;

% colors = [ % same strength
% 	153 185 255
% 	255 255 153
% 	153 255 168
% 	202 202 202
% ]/255;

colors = [
    210 255 255
    255 255 210
    255 210 255
    255 255 255
    ]/255;


%% make legends
ng = size(table.Pg,1);
genLabels = cell(1,ng);
genLabelsShort = cell(1,ng);
for i=1:ng
    strIdx = num2str(table2array(table.Pg(i,'GEN')));
    %genLabels{i} = ['GEN ' strIdx];
    busStr = num2str(table2array(table.Pg(i,'BUS')));
    genLabels{i} =  ['GEN ' strIdx ' BUS ' busStr];
    genLabelsShort{i} = num2str(table2array(table.Pg(i,'GEN')));
end

ncur = size(table.Beta,1);
curLabels = cell(1,ncur);
for i=1:ncur
    strIdx = num2str(table2array(table.Beta(i,'GEN')));
    busIdx = num2str(table2array(table.Beta(i,'BUS')));
    curLabels{i} = ['GEN ' strIdx ' BUS ' busIdx];
end

if plots.Pg
    %% Pg for all generators
    F1 = figure;
    cp = F1.Position;
    set(F1,'Position',[cp(1:2) fwidth fheight]);
    PgMat = table2array(table.Pg(:, setdiff(table.Pg.Properties.VariableNames,{'GEN','BUS','MIN','MAX'})));
    PH1 = bar(barx,PgMat');
    
    generator_greyscale(PH1,optns,table,color_factor,colors);
    
    grid on;
    %xlabel('Scenario');
    ylabel('P (MW)')
    columnlegend(n_legend_columns,genLabels,'Location','NorthWest');
    if plot_titles
        title(['Active power' '\newline' strrep(optns.caseName,'_',' ')])
    end
    yl = ylim;
    for i=1:length(PH1)
        tx = text([1:N]+PH1(i).XOffset-xoffset, PH1(i).YData+yoffset*(yl(2)-yl(1)), genLabelsShort{i},'FontSize',busFontSize);
        set(tx,'Rotation',90);
    end
    
    %ypos = 1800; % something below the maximum y-value
    ypos = yl(2);
    for i=1:4
        if i == 1 % start to the left
            xpos = 1+PH1(1).XOffset; % start at leftmost bar
            %xpos = 1+PH1(end).XOffset;
        else
            e = get(ltx,'Extent');
            %xpos = xpos + e(3); % left + width
            ypos = ypos - e(4) - get(ltx,'Margin'); % top - height
        end
        ltx = text(xpos,ypos,ann_strs{i},...
            'HorizontalAlignment','left','VerticalAlignment','top',...
            'BackgroundColor',color_factor(i)*colors(i,:),...
            'Margin',tx_margin);
    end
    
    for i=1:2
        e = get(ltx,'Extent');
        if i == 1
            ypos = ypos - 2*e(4);
        else
            ypos = ypos - e(4);
        end
        ltx = text(xpos,ypos,ann_strs2{i},...
            'HorizontalAlignment','left','VerticalAlignment','top',...
            'BackgroundColor',legend_greyfactor*ones(1,3),...
            'EdgeColor',[0 0 0],'LineStyle',edge_styles{i},...
            'LineWidth',1,'Margin',tx_margin);
    end
    set(F1,'Position',[cp(1:2) fwidth fheight+height_adjust]);
    set(gca,'FontSize',12);
    tight_axis(F1);
    %MaximizeFigureWindow();
    
    figureArray = [figureArray F1];
end


%% Qg for all generators
if plots.Qg
    F2 = figure;
    cp = F2.Position;
    set(F2,'Position',[cp(1:2) fwidth fheight]);
    PH2 = bar(barx,table2array(table.Qg(:, ...
        setdiff(table.Qg.Properties.VariableNames,{'GEN','BUS','MIN','MAX'})))' ...
        );
    
    generator_greyscale(PH2,optns,table,color_factor,colors);
    
    grid on;
    %xlabel('Scenario');
    ylabel('Q (MW)')
    %columnlegend(n_legend_columns,genLabels,'Location','NorthWest');
    if plot_titles
        title(['Reactive power' '\newline' strrep(optns.caseName,'_',' ')])
    end
    
    yl = ylim;
    for i=1:length(PH2)
        tx = text([1:N]+PH2(i).XOffset-xoffset, max(PH2(i).YData,zeros(1,N))+yoffset*(yl(2)-yl(1)), genLabelsShort{i},'FontSize',busFontSize);
        set(tx,'Rotation',90);
    end
    
    
    %ypos = 320; % something below the maximum y-value
    ypos = yl(2);
    for i=1:4
        if i == 1 % start to the left
            %xpos = 1+PH2(1).XOffset; % start at leftmost bar
            xpos = 1+PH2(end).XOffset;
        else
            e = get(ltx,'Extent');
            %xpos = xpos + e(3); % left + width
            ypos = ypos - e(4) - get(ltx,'Margin'); % top - height
        end
        ltx = text(xpos,ypos,ann_strs{i},...
            'HorizontalAlignment','left','VerticalAlignment','top',...
            'BackgroundColor',color_factor(i)*colors(i,:),...
            'Margin',tx_margin);
    end
    
    for i=1:2
        e = get(ltx,'Extent');
        if i == 1
            ypos = ypos - 2*e(4);
        else
            ypos = ypos - e(4);
        end
        ltx = text(xpos,ypos,ann_strs2{i},...
            'HorizontalAlignment','left','VerticalAlignment','top',...
            'BackgroundColor',legend_greyfactor*ones(1,3),...
            'EdgeColor',[0 0 0],'LineStyle',edge_styles{i},...
            'LineWidth',1,'Margin',tx_margin);
    end
    set(F2,'Position',[cp(1:2) fwidth fheight+height_adjust]);
    set(gca,'FontSize',12);
    tight_axis(F2);
    
    figureArray = [figureArray F2];
end

%% Voltages
if plots.Vm
    F3 = figure;
    cp = F3.Position;
    set(F3,'Position',[cp(1:2) fwidth fheight]);
    VmMat = table2array(table.Vm(:, setdiff(table.Vm.Properties.VariableNames,{'BUS','MIN','MAX'})));
    PH9 = bar(barx,VmMat');
    
    % categories:
    % No generation, Wind, Conventional, Wind+Conventional: 1,2,3,4
    catColors = [
        0 0 0; % black
        0 1 0; % green
        1 0 0; % red
        1 1 0 % yellow
        ];
    vmLabels = {};
    barLabels = {};
    for i=1:size(VmMat,1)
        vmLabels{i} = ['Bus ' num2str( table2array( table.Vm(i,'BUS') ))];
        barLabels{i} = num2str( table2array( table.Vm(i,'BUS') ) );
    end
    %legend(vmLabels);
    %columnlegend(n_legend_columns,vmLabels,'Location','NorthWest');
    
    voltage_greyscale(PH9,optns,table,color_factor,colors);
    
    
    grid on;
    %xlabel('Scenario');
    ylabel('V (pu)')
    if plot_titles
        title(['Voltage' '\newline' strrep(optns.caseName,'_',' ')]);
    end
    ylim([0.9 1.15]);
    yl = ylim;
    
    for i=1:length(PH9)
        tx = text([1:N]+PH9(i).XOffset-xoffset, PH9(i).YData+yoffset*(yl(2)-yl(1)), barLabels{i},'FontSize',busFontSize);
        set(tx,'Rotation',90);
    end
    
    rec_dim = [0.1 0.9 0.025 0.02];
    txt_dim = [0.13 0.85 0.1 0.1];
    
    
    
    %xpos = 1+PH9(1).XOffset; % start at leftmost bar
    xpos = 1+PH9(end).XOffset;
    for i=1:4
        if i == 1 % start to the left
            ypos = yl(2); % something below the maximum y-vale
        else
            e = get(ltx,'Extent');
            %xpos = xpos + e(3); % left + width
            ypos = ypos - e(4);
        end
        ltx = text(xpos,ypos,ann_strs{i},...
            'HorizontalAlignment','left','VerticalAlignment','top',...
            'BackgroundColor',color_factor(i)*colors(i,:),...
            'Margin',tx_margin);
    end
    
    set(F3,'Position',[cp(1:2) fwidth fheight+height_adjust]);
    set(gca,'FontSize',12);
    tight_axis(F3);
    
    figureArray = [figureArray F3];
end

%% Curtailment for Pcur
if plots.Curtail
    F4 = figure;
    PH4 = bar(1:N, table2array(table.Curtail(:, ...
        setdiff(table.Curtail.Properties.VariableNames,{'GEN','BUS'})))' ...
        );
    grid on;
    xlabel('Scenario');
    ylabel('Curtailed energy (MW)');
    legend(curLabels);
    if plot_titles
        title(['Curtailment' '\newline' strrep(optns.caseName,'_',' ')])
    end
    figureArray = [figureArray F4];
end
%% Beta for Pcur
if plots.Beta
    F5 = figure;
    PH5 = bar(1:N, table2array(table.Beta(:, ...
        setdiff(table.Beta.Properties.VariableNames,{'GEN','BUS'})))' ...
        );
    grid on;
    xlabel('Scenario');
    ylabel('Beta');
    legend(curLabels);
    if plot_titles
        title(['Curtailment fraction beta' '\newline' strrep(optns.caseName,'_',' ')])
    end
    figureArray = [figureArray F5];
end
%% Total curtailment

if plots.CurtailTotal
    curtail = table2array(table.Curtail(:,setdiff( table.Curtail.Properties.VariableNames,{'GEN','BUS'})));
    curtailedTotal = sum(curtail,1);
    
    F6 = figure;
    PH6 = bar(1:N, curtailedTotal);
    grid on;
    xlabel('Scenario');
    ylabel('Curtailment (MW)');
    if plot_titles
        title(['Total curtailment' '\newline' strrep(optns.caseName,'_',' ')])
    end
    figureArray = [figureArray F6];
end

%% Total curtailment in %
if plots.CurtailTotalPercent
    wind = table2array(table.Wind(:,setdiff(table.Wind.Properties.VariableNames,{'GEN','BUS'})));
    windTotal = sum(wind,1);
    
    F7 = figure;
    PH7 = bar(1:N, 100*curtailedTotal./windTotal);
    grid on;
    xlabel('Scenario');
    ylabel('Curtailment (%)');
    if plot_titles
        title(['Total curtailment' '\newline' strrep(optns.caseName,'_',' ')])
    end
    figureArray = [figureArray F7];
end
%% Plot transfers in corridors
if plots.Transfer
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
    if plot_titles
        title(['Transfer through corridors' '\newline' strrep(optns.caseName,'_',' ')]);
    end
    figureArray = [figureArray F8];
    
end

%% Wind Penetration
if plots.Penetration
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
    
    F9 = figure;
    PH9 = bar(1:N,penetration*100);
    grid on;
    xlabel('Scenario');
    ylabel('Wind Penetration (%)');
    if plot_titles
        title(['Penetration' '\newline' strrep(optns.caseName,'_',' ')]);
    end
    figureArray = [figureArray F9];
end

%% save figures

% figureArray = [F1 F2 F3 F4 F5 F6 F7 F8 F3];
% figNames = {'Pg','Qg','WindPenetration','CutailMW','CurtailFraction','CurtailTot','CurtailTotPc','Transfers','Vm'};
% if optns.saveFigures
%    for i=1:length(figureArray)
%        thisFigureName = [optns.caseName '_' figNames{i}];
%         saveas(figureArray(i),['figures/' thisFigureName '.fig']);
%         saveas(figureArray(i),['figures/' thisFigureName '.png']);
%    end
% end

