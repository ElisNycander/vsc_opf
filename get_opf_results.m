function [results,tab] = get_opf_results(om,x,Lambda,optns)
% Collect results from opf into matrices with one column for each
% contingency.
% Concerning indexing: 
% Internal indexing for buses changes BUS_I but not the order of the rows
% Internal indexing for generators only changes order of the rows

%% for script
% clear;
% close all;
% load('vscopf_post_processing_data.mat');



%% Notes
% 1. To change the table columns, they must be changed in this file, in 
% the following places:
% * Making the table data
% * Making the table labels
% The column indices must also be updated in 
% * plot_results
% * verify_solutions
%
%%
define_constants;
[vv, ll, nn] = get_idx(om);
[v0,vl,vu] = getv(om);
mpc = get_mpc(om);
mpopt = optns.mpopt;
[cs,gen,bus,bus2,gen2,branch,baseMVA] = deal(mpc.contingencies,mpc.gen,mpc.bus,mpc.bus2,mpc.gen2,mpc.branch,mpc.baseMVA);

nc = cs.N;
ng = size(gen,1);
nb = size(bus,1);
nl = size(branch,1);

[~,Yf,Yt] = makeYbus(baseMVA,bus,branch);

numIdxCurtail = find(gen2(:,PTYPE)==PCUR); % numbered index
idxNotCurtail = find(gen2(:,PTYPE)~=PCUR);
%bolIdxCurtail = gen2(:,PTYPE) == PCUR;
nCurtail = length(numIdxCurtail);
%% PRINT GENERATION

% collect injections in matrix:
%     P0  P1  ... Pnc (cs)
% G1
% G2
% ...
% Gng 
% (generators)


%% containers
Pg = zeros(ng,nc);
Qg = zeros(ng,nc);
% lagrange multipliers
PgU = zeros(ng,nc);
PgL = zeros(ng,nc);
QgU = zeros(ng,nc);
QgL = zeros(ng,nc);

Vm = zeros(nb,nc);
Va = zeros(size(Vm));
VmU = zeros(size(Vm));
VmL = zeros(size(Vm));
VaU = zeros(size(Vm));
VaL = zeros(size(Vm));

Beta = zeros(nCurtail,nc);
Curtail = zeros(nCurtail,nc);
Wind = zeros(nCurtail,nc);

flow = nan(nl,2*nc);
Sflow = nan(nl,2*nc);
Pflow = nan(nl,2*nc);
Iflow = nan(nl,2*nc);
Qflow = nan(nl,2*nc);
Ploss = nan(nl,2*nc);

flow_lam = nan(size(flow));
% Pg_lamU = nan(size(Pg));
% Pg_lamL = nan(size(Pg));
% Vm_lamU = nan(size(Vm));
% Vm_lamL = nan(size(Vm));

corrFrom = zeros(length(optns.transferCorridors),nc); 
corrTo = zeros(length(optns.transferCorridors),nc);


countLines = 0;
countConstrainedLines = 1;

includePgBase = isfield(vv.i1,'Pg');
if includePgBase
    fixBaseWind = vv.N.Pg < ng;
end
% if isfield(vv.i1,'Pg')
%     includePgBase = 1;
% else
%     includePgBase = 0;
% end
for i=1:nc
    % indices 
    if i == 1
        si = '';
    else 
        si = num2str(i);
    end
    
    iVm = vv.i1.(['Vm' si]):vv.iN.(['Vm' si]);
    iVa = vv.i1.(['Va' si]):vv.iN.(['Va' si]);
    if i > 1 || includePgBase
        iPg = vv.i1.(['Pg' si]):vv.iN.(['Pg' si]);
    else
        iPg = [];
    end
    iQg = vv.i1.(['Qg' si]):vv.iN.(['Qg' si]);
    if i > 1
        iBeta = vv.i1.(['Beta' si]):vv.iN.(['Beta' si]);
    end
    % voltages
    Vm(:,i) = x(iVm);
    Va(:,i) = x(iVa);
    V = x(iVm) .* exp(1j * x(iVa));
    % lagrange multipliers
    %Vm_lamU(:,i) = Lambda.upper((['Vm' si]):vv.iN.(['Vm' si]));
    %Vm_lamL(:,i) = Lambda.lower((['Vm' si]):vv.iN.(['Vm' si]));
    
    % generation
    idxTrip = ~cs.activeGenerators(:,i);
    
    %idxPfix = mpc.gen2(:,PTYPE) == PFIX;
    idxPvar = and(mpc.gen2(:,PTYPE) == PVAR,~idxTrip);
    idxPfix = mpc.gen2(:,PTYPE) == PFIX;
    idxPcur = mpc.gen2(:,PTYPE) == PCUR;
    idxQvar = ~idxTrip;
    
    
    %idxPvar = and(~idxTrip,idxPVar);
    %idxQvar = and(~idxTrip,idxQVar);
%     idxActiveGenerators = cs.activeGenerators(:,i);
%     idxTrippedGenerators = ~cs.activeGenerators(:,i);
%     nActiveGenerators = cs.nActiveGenerators(i);
%     
    if i > 1
        Pg(idxPvar,i) = x(iPg);
        Pg(idxPfix,i) = Pg(idxPfix,1); % store fixed values in all contingencies
        Pg(idxPcur,i) = cs.wind(1+(i-1)*nCurtail:i*nCurtail)/baseMVA.*(1-x(iBeta)); % store wind scenarios
        Pg(idxTrip,i) = NaN;
        PgL(idxPvar,i) = Lambda.lower(iPg);
        PgU(idxPvar,i) = Lambda.upper(iPg);
        Qg(idxQvar,i) = x(iQg);
        Qg(~idxQvar,i) = Qg(~idxQvar,1); % store fixed values in all contingencies
        Qg(idxTrip,i) = NaN;
        QgL(idxQvar,i) = Lambda.lower(iQg);
        QgU(idxQvar,i) = Lambda.upper(iQg);
        Beta(:,i) = x(iBeta);
        Curtail(:,i) = cs.wind(1+(i-1)*nCurtail:i*nCurtail).*x(iBeta);
        
    else % base case
        if includePgBase
            if ~fixBaseWind % all Pg are optimization variables
                Pg(:,i) = x(iPg);
                PgL(:,i) = Lambda.lower(iPg);
                PgU(:,i) = Lambda.upper(iPg);
            else % only non-curtailable Pg are optimization variables
                Pg(idxNotCurtail,i) = x(iPg);
                PgL(idxNotCurtail,i) = Lambda.lower(iPg);
                PgU(idxNotCurtail,i) = Lambda.upper(iPg);
                
                % get pre-set Pg values for curtailable generators
                Pg(numIdxCurtail,i) = mpc.gen(numIdxCurtail,PG)/mpc.baseMVA;
                PgL(numIdxCurtail,i) = NaN;
                PgU(numIdxCurtail,i) = NaN;
            end
        else
            Pg(:,i) = mpc.gen(:,PG)/mpc.baseMVA; % get pre-set Pg values
        end
        Beta(:,i) = 0;
        Curtail(:,i) = 0;
        Qg(:,i) = x(iQg);
        QgL(:,i) = Lambda.lower(iQg);
        QgU(:,i) = Lambda.upper(iQg);
        
    end
    
    Wind(:,i) = cs.wind(1+(i-1)*nCurtail:i*nCurtail); 
    % voltage limits
    VmU(:,i) = Lambda.upper(iVm);
    VmL(:,i) = Lambda.lower(iVm);
    VaU(:,i) = Lambda.upper(iVa);
    VaL(:,i) = Lambda.lower(iVa);
    

    % calculate branch flows
    active_lines = cs.activeLines(:,i);
    inl = cs.nActiveLines(i);
    
    iYf = cs.Yf(countLines+1:countLines+inl,:);
    iYt = cs.Yt(countLines+1:countLines+inl,:);
    
    
    If = iYf * V;
    It = iYt * V;
    Iflow(active_lines,2*i-1) = If.*conj(If);
    Iflow(active_lines,2*i) = It.*conj(It);
    
    iSf = V(branch(active_lines, F_BUS)) .* conj(iYf * V);  %% complex power injected at "from" bus (p.u.)
    iSt = V(branch(active_lines, T_BUS)) .* conj(iYt * V);  %% complex power injected at "to" bus (p.u.)
    
    
    Pflow(active_lines,2*i-1) = real(iSf)*baseMVA;
    Pflow(active_lines,2*i) = real(iSt)*baseMVA;
    Ploss(active_lines,2*i-1) = ( real(iSf) + real(iSt) ) * baseMVA;
    
    Pinj = zeros(size(iSf));
    for ii=1:length(iSf) 
        Pinj(ii) = max( real(iSf(ii)),real(iSt(ii)) ); 
    end
    Ploss(active_lines,2*i) = Ploss(active_lines,2*i-1) ./ Pinj;
    
    Qflow(active_lines,2*i-1) = imag(iSf)*baseMVA;
    Qflow(active_lines,2*i) = imag(iSt)*baseMVA;
    
    Sflow(active_lines,2*i-1) = (abs(iSf))*baseMVA;
    Sflow(active_lines,2*i) = (abs(iSt))*baseMVA;
    


    % lagrange multipliers for line constraints
    if optns.gen.optimizeBaseP
    flow_lam(cs.constrainedActiveLines(:,i),2*i-1) = Lambda.ineqnonlin(...
            countConstrainedLines : countConstrainedLines + cs.nConstrainedActiveLines(i) - 1 );
    flow_lam(cs.constrainedActiveLines(:,i),2*i) = Lambda.ineqnonlin(...
            countConstrainedLines+cs.nConstrainedActiveLines(i) : ...
            countConstrainedLines+2*cs.nConstrainedActiveLines(i) - 1);
    else
        if i == 1
            flow_lam(cs.constrainedActiveLines(:,i),2*i-1) = NaN;
            flow_lam(cs.constrainedActiveLines(:,i),2*i) = NaN;
        else
            flow_lam(cs.constrainedActiveLines(:,i),2*i-1) = Lambda.ineqnonlin(...
                countConstrainedLines : countConstrainedLines + cs.nConstrainedActiveLines(i) - 1 );
            flow_lam(cs.constrainedActiveLines(:,i),2*i) = Lambda.ineqnonlin(...
                countConstrainedLines+cs.nConstrainedActiveLines(i) : ...
                countConstrainedLines+2*cs.nConstrainedActiveLines(i) - 1);
        end
        
        
    end
    
    % sum flows over corridors
    for j=1:length(optns.transferCorridors)
        thisCorridor = optns.corridorIdxs{j};
        direction = sign(thisCorridor);
        thisCorridor = abs(thisCorridor);
        
        for jj=1:length(thisCorridor)
            if ~isnan(Pflow(thisCorridor(jj),2*i))
                if direction(jj) > 0
                    corrFrom(j,i) = corrFrom(j,i) + Pflow(thisCorridor(jj),2*i-1);
                    corrTo(j,i) = corrTo(j,i) + Pflow(thisCorridor(jj),2*i);
                else % switch from and to bus for this line
                    corrFrom(j,i) = corrFrom(j,i) + Pflow(thisCorridor(jj),2*i);
                    corrTo(j,i) = corrTo(j,i) + Pflow(thisCorridor(jj),2*i-1);
                end
            end
        end
    end
    
    %% increment counters
    if optns.gen.optimizeBaseP || i > 1
        countConstrainedLines = countConstrainedLines + 2*cs.nConstrainedActiveLines(i);
    end
    countLines = countLines + inl;
end

% choose which type of flow are used for ratings
if strcmp(mpopt.opf.flow_lim,'I')
    flow = Iflow;
elseif strcmp(mpopt.opf.flow_lim,'P')
    flow = Pflow;
else
    flow = Sflow;
end




% put gen bus in matrix, convert to nominal units
Pg = [mpc.order.gen.e2i mpc.order.bus.i2e(mpc.gen(:,GEN_BUS)) mpc.gen(:,PMIN) mpc.gen(:,PMAX) Pg*mpc.baseMVA];
Qg = [mpc.order.gen.e2i mpc.order.bus.i2e(mpc.gen(:,GEN_BUS)) mpc.gen(:,QMIN) mpc.gen(:,QMAX) Qg*mpc.baseMVA];

genExt = mpc.order.gen.e2i;
busExt = mpc.order.bus.i2e(mpc.gen(:,GEN_BUS));


PgU = [genExt busExt PgU];
PgL = [genExt busExt PgL];
QgU = [genExt busExt QgU];
QgL = [genExt busExt QgL];

%Beta = [mpc.order.gen.i2e(idxCurtail) mpc.order.bus.i2e(mpc.gen(idxCurtail,GEN_BUS)) Beta];
Beta = [mpc.order.gen.e2i(numIdxCurtail) mpc.order.bus.i2e(mpc.gen(numIdxCurtail,GEN_BUS)) Beta];
expCurtail = [cs.probabilities(1:end)'; sum(Curtail,1); ...
              cs.probabilities(1:end)'.*sum(Curtail,1)];
Curtail = [mpc.order.gen.e2i(numIdxCurtail) mpc.order.bus.i2e(mpc.gen(numIdxCurtail,GEN_BUS)) Curtail];
Wind = [mpc.order.gen.e2i(numIdxCurtail) mpc.order.bus.i2e(mpc.gen(numIdxCurtail,GEN_BUS)) Wind];


s = mpc.order.gen.i2e;
% sort according to external indexing
Pg = Pg(s,:);
Qg = Qg(s,:);
PgU = PgU(s,:);
PgL = PgL(s,:);
QgU = QgU(s,:);
QgL = QgL(s,:);

[~,sortIdx] = sort(mpc.order.gen.e2i(numIdxCurtail));
Beta = Beta(sortIdx,:);
Curtail = Curtail(sortIdx,:);
Wind = Wind(sortIdx,:);
% 

% put bus nr in Va
Va = [mpc.order.bus.i2e(mpc.bus(:,BUS_I)) Va*180/pi];
Vm = [mpc.order.bus.i2e(mpc.bus(:,BUS_I)) mpc.bus(:,VMIN) mpc.bus(:,VMAX) Vm];
busExt = mpc.order.bus.i2e(mpc.bus(:,BUS_I));

VmU = [busExt VmU];
VmL = [busExt VmL];
VaU = [busExt VaU];
VaL = [busExt VaL];

flow = [(1:size(branch,1))' mpc.order.bus.i2e(branch(:,F_BUS)) mpc.order.bus.i2e(branch(:,T_BUS)) branch(:,RATE_A) flow];
flow_lam = [(1:size(branch,1))' mpc.order.bus.i2e(branch(:,F_BUS)) mpc.order.bus.i2e(branch(:,T_BUS)) branch(:,RATE_A) flow_lam];


Pflow = [(1:size(branch,1))' mpc.order.bus.i2e(branch(:,F_BUS)) mpc.order.bus.i2e(branch(:,T_BUS)) Pflow];
Qflow = [(1:size(branch,1))' mpc.order.bus.i2e(branch(:,F_BUS)) mpc.order.bus.i2e(branch(:,T_BUS)) Qflow];
Ploss = [(1:size(branch,1))' mpc.order.bus.i2e(branch(:,F_BUS)) mpc.order.bus.i2e(branch(:,T_BUS)) Ploss];

corrFrom = [(1:size(corrFrom,1))' corrFrom];
corrTo = [(1:size(corrTo,1))' corrTo];

% set multipliers below threshold to 0
thrs = optns.lambdaTolerance;
flow_lam(flow_lam < thrs) = 0;
PgU(PgU < thrs) = 0;
PgL(PgL < thrs) = 0;
QgU(QgU < thrs) = 0;
QgL(QgL < thrs) = 0;
VmU(VmU < thrs) = 0;
VmL(VmL < thrs) = 0;
VaU(VaU < thrs) = 0;
VaL(VaL < thrs) = 0;

%% make tables
varnames_beta = {'GEN', 'BUS'};
varnames_curtail = {'GEN','BUS'};
varnames_pg = {'GEN', 'BUS' 'MIN' 'MAX'};
varnames_qg = {'GEN', 'BUS' 'MIN' 'MAX'};
varnames_pglam = {'GEN','BUS'};
varnames_qglam = {'GEN','BUS'};
varnames_va = {'BUS'};
varnames_vm = {'BUS','MIN','MAX'};
varnames_S = {'BRANCH','FROM','TO','RATE_A'};
varnames_expcurtail = {};
varnames_wind = {'GEN','BUS'};
varnames_vmlam = {'BUS'};
varnames_corridors = {'CORRIDOR'};
varnames_pflow = {'BRANCH','FROM','TO'};
varnames_qflow = {'BRANCH','FROM','TO'};
varnames_ploss = {'BRANCH','FROM','TO'};
for i=1:nc
	if i == 1
		stringIdx = '';
	else
		stringIdx = num2str(i);
	end
	
    varnames_pg{i+4} = ['PG' stringIdx];
    varnames_qg{i+4} = ['QG' stringIdx];
    varnames_pglam{i+2} = ['PG' stringIdx];
    varnames_qglam{i+2} = ['QG' stringIdx];
    varnames_va{i+1} = ['VA' stringIdx];
    varnames_vmlam{i+1} = ['VM' stringIdx];
    varnames_vm{i+3} = ['VM' stringIdx];
    varnames_S{2*i+3} = ['SF' stringIdx];
    varnames_S{2*i+4} = ['ST' stringIdx];
    varnames_beta{i+2} = ['BETA' stringIdx];
    varnames_curtail{i+2} = ['PGC' stringIdx];
    varnames_wind{i+2} = ['PWIND' stringIdx];
    varnames_expcurtail{i} = ['EPGC' stringIdx];
    varnames_corridors{i+1} = ['T' stringIdx];
    varnames_pflow{2*i+2} = ['PF' stringIdx];
    varnames_pflow{2*i+3} = ['PT' stringIdx];
    varnames_qflow{2*i+2} = ['QF' stringIdx];
    varnames_qflow{2*i+3} = ['QT' stringIdx];
    varnames_ploss{2*i+2} = ['MW' stringIdx];
    varnames_ploss{2*i+3} = ['PC' stringIdx];
end

Pg_table = array2table(Pg,...
    'VariableNames',varnames_pg);
Qg_table = array2table(Qg,...
    'VariableNames',varnames_qg);
Vm_table = array2table(Vm,'VariableNames',varnames_vm);
Va_table = array2table(Va,'VariableNames',varnames_va);


Pg_lamU_table = array2table(PgU,'VariableNames',varnames_pglam);
Pg_lamL_table = array2table(PgL,'VariableNames',varnames_pglam);
Qg_lamU_table = array2table(QgU,'VariableNames',varnames_qglam);
Qg_lamL_table = array2table(QgL,'VariableNames',varnames_qglam);
Vm_lamU_table = array2table(VmU,'VariableNames',varnames_vmlam);
Vm_lamL_table = array2table(VmL,'VariableNames',varnames_vmlam);
Va_lamU_table = array2table(VaU,'VariableNames',varnames_va);
Va_lamL_table = array2table(VaL,'VariableNames',varnames_va);

Beta_table = array2table(Beta,'VariableNames',varnames_beta);
Curtail_table = array2table(Curtail,'VariableNames',varnames_curtail);
expCurtail_table = array2table(expCurtail,'VariableNames',varnames_expcurtail);
Wind_table = array2table(Wind,'VariableNames',varnames_wind);
S_table = array2table(flow,'VariableNames',varnames_S);
S_lam_table = array2table(flow_lam,'VariableNames',varnames_S);
transferFrom_table = array2table(corrFrom,'VariableNames',varnames_corridors);
transferTo_table = array2table(corrTo,'VariableNames',varnames_corridors);

Pflow_table = array2table(Pflow,'VariableNames',varnames_pflow);
Qflow_table = array2table(Qflow,'VariableNames',varnames_qflow);
Ploss_table = array2table(Ploss,'VariableNames',varnames_ploss);

%% table with PQ-capability constraints
pqGens = find(optns.gen.pqFactor);
qCols = {};
pCols = {};
for i=1:nc
    if i == 1
        stringIdx = '';
    else
        stringIdx = num2str(i);
    end
    qCols{i}=['QG' stringIdx];
    pCols{i}=['PG' stringIdx];
    
    
end
Ptmp = Pg_table(pqGens,pCols);
Qtmp = Qg_table(pqGens,qCols);
Label_tmp = Pg_table(pqGens,{'GEN','BUS'});
Cap_tmp = array2table(optns.gen.pqFactor(pqGens),'VariableNames',{'PQ_FACTOR'});
PQ_table = [Label_tmp Cap_tmp Ptmp Qtmp];

%% find non-zero multipliers
str = '';
vars = {'PgU','PgL','QgU','QgL','VmU','VmL','VaU','VaL','flow_lam'};

for i=1:length(vars)
    switch vars{i} % exclude columns which don't contain multipliers
        case 'PgU'
            exclColIdx = [1 2];
        case 'PgL'
            exclColIdx = [1 2];
        case 'QgU'
            exclColIdx = [1 2];
        case 'QgL'
            exclColIdx = [1 2];
        case 'VmU'
            exclColIdx = [1];
        case 'VmL'
            exclColIdx = [1];
        case 'VaU'
            exclColIdx = [1];
        case 'VaL'
            exclColIdx = [1];
        case 'flow_lam'
            exclColIdx = 1:4;
    end
        
    eval(['mat=' vars{i} ';']);
    mat(:,exclColIdx) = []; % remove columns
    nzIdx = find(mat > optns.lambdaTolerance);
    nanIdx = find(isnan(mat));
    nzIdx = setdiff(nzIdx,nanIdx); % non-zero entries that are not NaN
    
    if ~isempty(nzIdx)
        % add information to string
        for j=1:length(nzIdx)
            [row,col] = ind2sub(size(mat),nzIdx(j));
            sInfo = '';
            switch vars{i}
                case 'PgU'
                    sInfo = sprintf('PG max: Gen %i, Scenario %i',[PgU(row,1) col]);
                case 'PgL'
                    sInfo = sprintf('PG min: Gen %i, Scenario %i',[PgL(row,1) col]);
                case 'QgU'
                     sInfo = sprintf('QG max: Gen %i, Scenario %i',[QgU(row,1) col]);
                case 'QgL'
                    sInfo = sprintf('QG min: Gen %i, Scenario %i',[QgL(row,1) col]);
                case 'VmU'
                   sInfo = sprintf('VM max: Bus %i, Scenario %i',[VmU(row,1) col]);
                case 'VmL'
                    sInfo = sprintf('VM min: Bus %i, Scenario %i',[VmL(row,1) col]);
                case 'VaU'
                    sInfo = sprintf('VA max: Bus %i, Scenario %i',[VaU(row,1) col]);
                case 'VaL'
                    sInfo = sprintf('VA min: Bus %i, Scenario %i',[VaL(row,1) col]);
                case 'flow_lam'
                    sInfo = sprintf('flow Max: Branch %i, %i-%i, Scenario %i',[flow_lam(row,[1:3]) ceil(col/2)]);
                    
            end
            str = [str sprintf('\n') sInfo];
        end
    end
end
lamInfo = str;

tab = struct();
[tab.Pg,tab.Qg,tab.Va,tab.Vm, tab.Flow, tab.Flowlam, tab.Curtail] = deal(Pg_table,Qg_table,Va_table,Vm_table, S_table, S_lam_table, Curtail_table);
[tab.PgUlam,tab.PgLlam,tab.QgUlam,tab.QgLlam,tab.VmUlam,tab.VmLlam,tab.VaUlam,tab.VaLlam, tab.Beta] = ...
    deal(Pg_lamU_table,Pg_lamL_table,Qg_lamU_table,Qg_lamL_table,Vm_lamU_table,Vm_lamL_table,Va_lamU_table,Va_lamL_table, Beta_table);
[tab.ExpCurtail, tab.Wind, tab.PQ] = deal(expCurtail_table, Wind_table,PQ_table);
[tab.lamInfo, tab.transferFrom, tab.transferTo, tab.Pflow, tab.Qflow] = deal(lamInfo, transferFrom_table, transferTo_table, Pflow_table, Qflow_table);
[tab.Ploss] = deal(Ploss_table);
results = mpc;
[results.bus, results.branch, results.gen, results.gen2, results.bus2, ...
    results.om] = ...
        deal(bus, branch, gen, gen2, bus2, om);
    
    
    % calculate branch flows through corridors


% for i=1:size(corrFrom,1)
%     % find indices of lines for this corridor
%     thisCorridor = optns.corridorIdxs{i};
%     
%     direction = sign(thisCorridor);
%     thisCorridor = abs(thisCorridor);
%     
%     for j=1:nc
%   
%         
%     end
% end
    

