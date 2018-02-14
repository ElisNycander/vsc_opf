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

%%
define_constants;
[vv, ll, nn] = get_idx(om);
mpc = get_mpc(om);
mpopt = optns.mpopt;
[cs,gen,bus,bus2,gen2,branch,baseMVA] = deal(mpc.contingencies,mpc.gen,mpc.bus,mpc.bus2,mpc.gen2,mpc.branch,mpc.baseMVA);

nc = cs.N;
ng = size(gen,1);
nb = size(bus,1);
nl = size(branch,1);

[~,Yf,Yt] = makeYbus(baseMVA,bus,branch);

idxCurtail = find(gen2(:,PTYPE)==PCUR);
nCurtail = length(idxCurtail);
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

Sflow = nan(nl,2*nc);

Sflow_lam = nan(size(Sflow));
Pg_lamU = nan(size(Pg));
Pg_lamL = nan(size(Pg));
Vm_lamU = nan(size(Vm));
Vm_lamL = nan(size(Vm));

countLines = 0;
countConstrainedLines = 1;

if isfield(vv.i1,'Pg')
    includePgBase = 1;
else
    includePgBase = 0;
end
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
    idxActiveGenerators = cs.activeGenerators(:,i);
    idxTrippedGenerators = ~cs.activeGenerators(:,i);
    nActiveGenerators = cs.nActiveGenerators(i);
    
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
            Pg(:,i) = x(iPg);
            PgL(:,i) = Lambda.lower(iPg);
            PgU(:,i) = Lambda.upper(iPg);
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
    
    
    if strcmp(mpopt.opf.flow_lim,'I')
        If = iYf * V;
        It = iYt * V;
        Sflow(active_lines,2*i-1) = If.*conj(If);
        Sflow(active_lines,2*i) = It.*conj(It);
    else
        iSf = V(branch(active_lines, F_BUS)) .* conj(iYf * V);  %% complex power injected at "from" bus (p.u.)
        iSt = V(branch(active_lines, T_BUS)) .* conj(iYt * V);  %% complex power injected at "to" bus (p.u.)
        
        if strcmp(mpopt.opf.flow_lim,'P')
            %Sf(:,i) = real(iSf)*baseMVA;
            Sflow(active_lines,2*i-1) = real(iSf)*baseMVA;
            Sflow(active_lines,2*i) = real(iSt)*baseMVA;
        else
            %Sf(:,i) = abs(iSf)*baseMVA;
            Sflow(active_lines,2*i-1) = (abs(iSf)).^2*baseMVA;
            Sflow(active_lines,2*i) = (abs(iSt)).^2*baseMVA;
        end
        
    end
    % lagrange multipliers for line constraints
    Sflow_lam(cs.constrainedActiveLines(:,i),2*i-1) = Lambda.ineqnonlin(...
            countConstrainedLines : countConstrainedLines + cs.nConstrainedActiveLines(i) - 1 );
    Sflow_lam(cs.constrainedActiveLines(:,i),2*i) = Lambda.ineqnonlin(...
            countConstrainedLines+cs.nConstrainedActiveLines(i) : ...
            countConstrainedLines+2*cs.nConstrainedActiveLines(i) - 1);
        
    %% increment counters    
    countConstrainedLines = countConstrainedLines + 2*cs.nConstrainedActiveLines(i);
    countLines = countLines + inl;
end

% put gen bus in matrix, convert to nominal units
Pg = [mpc.order.gen.e2i mpc.order.bus.i2e(mpc.gen(:,GEN_BUS)) Pg*mpc.baseMVA];
Qg = [mpc.order.gen.e2i mpc.order.bus.i2e(mpc.gen(:,GEN_BUS)) Qg*mpc.baseMVA];

genExt = mpc.order.gen.e2i;
busExt = mpc.order.bus.i2e(mpc.gen(:,GEN_BUS));

PgU = [genExt busExt PgU];
PgL = [genExt busExt PgL];
QgU = [genExt busExt QgU];
QgL = [genExt busExt QgL];

%Beta = [mpc.order.gen.i2e(idxCurtail) mpc.order.bus.i2e(mpc.gen(idxCurtail,GEN_BUS)) Beta];
Beta = [mpc.order.gen.e2i(idxCurtail) mpc.order.bus.i2e(mpc.gen(idxCurtail,GEN_BUS)) Beta];
expCurtail = [cs.probabilities(1:end)'; sum(Curtail,1); ...
              cs.probabilities(1:end)'.*sum(Curtail,1)];
Curtail = [mpc.order.gen.e2i(idxCurtail) mpc.order.bus.i2e(mpc.gen(idxCurtail,GEN_BUS)) Curtail];
Wind = [mpc.order.gen.e2i(idxCurtail) mpc.order.bus.i2e(mpc.gen(idxCurtail,GEN_BUS)) Wind];
% add row with sum of curtailment
% Curtail = [Curtail; sum(Curtail,1)]; 
% Curtail(end,1:2) = NaN;
% table with sum of curtailment and probabilities
%expCurtail = sum(Curtail(:,3:end),1)


s = mpc.order.gen.i2e;
% sort according to external indexing
Pg = Pg(s,:);
Qg = Qg(s,:);
PgU = PgU(s,:);
PgL = PgL(s,:);
QgU = QgU(s,:);
QgL = QgL(s,:);

[~,sortIdx] = sort(mpc.order.gen.e2i(idxCurtail));
Beta = Beta(sortIdx,:);
Curtail = Curtail(sortIdx,:);
Wind = Wind(sortIdx,:);
% 

% put bus nr in Va
Va = [mpc.order.bus.i2e(mpc.bus(:,BUS_I)) Va*180/pi];
Vm = [mpc.order.bus.i2e(mpc.bus(:,BUS_I)) Vm];
busExt = mpc.order.bus.i2e(mpc.bus(:,BUS_I));

VmU = [busExt VmU];
VmL = [busExt VmL];
VaU = [busExt VaU];
VaL = [busExt VaL];

Sflow = [(1:size(branch,1))' mpc.order.bus.i2e(branch(:,F_BUS)) mpc.order.bus.i2e(branch(:,T_BUS)) branch(:,RATE_A) Sflow];
Sflow_lam = [(1:size(branch,1))' mpc.order.bus.i2e(branch(:,F_BUS)) mpc.order.bus.i2e(branch(:,T_BUS)) branch(:,RATE_A) Sflow_lam];

% set multipliers below threshold to 0
thrs = optns.lambdaTolerance;
Sflow_lam(Sflow_lam < thrs) = 0;
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
varnames_pg = {'GEN', 'BUS'};
varnames_qg = {'GEN', 'BUS'};
varnames_va = {'BUS'};
varnames_vm = {'BUS'};
varnames_S = {'BRANCH','FROM','TO','RATE_A'};
varnames_expcurtail = {};
varnames_wind = {'GEN','BUS'};
for i=1:nc
	if i == 1
		stringIdx = '';
	else
		stringIdx = num2str(i);
	end
	
    varnames_pg{i+2} = ['PG' stringIdx];
    varnames_qg{i+2} = ['QG' stringIdx];
    varnames_va{i+1} = ['VA' stringIdx];
    varnames_vm{i+1} = ['VM' stringIdx];
    varnames_S{2*i+3} = ['SF' stringIdx];
    varnames_S{2*i+4} = ['ST' stringIdx];
    varnames_beta{i+2} = ['BETA' stringIdx];
    varnames_curtail{i+2} = ['PGC' stringIdx];
    varnames_wind{i+2} = ['PWIND' stringIdx];
    varnames_expcurtail{i} = ['EPGC' stringIdx];
    
end

Pg_table = array2table(Pg,...
    'VariableNames',varnames_pg);
Qg_table = array2table(Qg,...
    'VariableNames',varnames_qg);
Vm_table = array2table(Vm,'VariableNames',varnames_vm);
Va_table = array2table(Va,'VariableNames',varnames_va);

Pg_lamU_table = array2table(PgU,'VariableNames',varnames_pg);
Pg_lamL_table = array2table(PgL,'VariableNames',varnames_pg);
Qg_lamU_table = array2table(QgU,'VariableNames',varnames_qg);
Qg_lamL_table = array2table(QgL,'VariableNames',varnames_qg);
Vm_lamU_table = array2table(VmU,'VariableNames',varnames_vm);
Vm_lamL_table = array2table(VmL,'VariableNames',varnames_vm);
Va_lamU_table = array2table(VaU,'VariableNames',varnames_va);
Va_lamL_table = array2table(VaL,'VariableNames',varnames_va);

Beta_table = array2table(Beta,'VariableNames',varnames_beta);
Curtail_table = array2table(Curtail,'VariableNames',varnames_curtail);
expCurtail_table = array2table(expCurtail,'VariableNames',varnames_expcurtail);
Wind_table = array2table(Wind,'VariableNames',varnames_wind);
S_table = array2table(Sflow,'VariableNames',varnames_S);
S_lam_table = array2table(Sflow_lam,'VariableNames',varnames_S);



%% find non-zero multipliers
str = '';
vars = {'PgU','PgL','QgU','QgL','VmU','VmL','VaU','VaL','Sflow_lam'};

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
        case 'Sflow_lam'
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
                case 'Sflow_lam'
                    sInfo = sprintf('Sflow Max: Branch %i, %i-%i, Scenario %i',[Sflow_lam(row,[1:3]) ceil(col/2)]);
                    
            end
            str = [str sprintf('\n') sInfo];
        end
    end
end
lamInfo = str;

tab = struct();
[tab.Pg,tab.Qg,tab.Va,tab.Vm, tab.S, tab.Slam, tab.Curtail] = deal(Pg_table,Qg_table,Va_table,Vm_table, S_table, S_lam_table, Curtail_table);
[tab.PgUlam,tab.PgLlam,tab.QgUlam,tab.QgLlam,tab.VmUlam,tab.VmLlam,tab.VaUlam,tab.VaLlam, tab.Beta] = ...
    deal(Pg_lamU_table,Pg_lamL_table,Qg_lamU_table,Qg_lamL_table,Vm_lamU_table,Vm_lamL_table,Va_lamU_table,Va_lamL_table, Beta_table);
[tab.ExpCurtail, tab.Wind] = deal(expCurtail_table, Wind_table);
[tab.lamInfo] = deal(lamInfo);

results = mpc;
[results.bus, results.branch, results.gen, results.gen2, results.bus2, ...
    results.om] = ...
        deal(bus, branch, gen, gen2, bus2, om);
