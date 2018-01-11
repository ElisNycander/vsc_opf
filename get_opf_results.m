function [results,tab] = get_opf_results(om,x,Lambda,optns)

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

Sflow = nan(nl,2*nc);

Sflow_lam = nan(size(Sflow));
Pg_lamU = nan(size(Pg));
Pg_lamL = nan(size(Pg));
Vm_lamU = nan(size(Vm));
Vm_lamL = nan(size(Vm));

countLines = 0;
countConstrainedLines = 1;
for i=1:nc
    % indices 
    if i == 1
        si = '';
    else 
        si = num2str(i);
    end
    iVm = vv.i1.(['Vm' si]):vv.iN.(['Vm' si]);
    iVa = vv.i1.(['Va' si]):vv.iN.(['Va' si]);
    iPg = vv.i1.(['Pg' si]):vv.iN.(['Pg' si]);
    iQg = vv.i1.(['Qg' si]):vv.iN.(['Qg' si]);
    
    % voltages
    Vm(:,i) = x(iVm);
    Va(:,i) = x(iVa);
    V = x(iVm) .* exp(1j * x(iVa));
    % lagrange multipliers
    %Vm_lamU(:,i) = Lambda.upper((['Vm' si]):vv.iN.(['Vm' si]));
    %Vm_lamL(:,i) = Lambda.lower((['Vm' si]):vv.iN.(['Vm' si]));
    
    % generation
    idxPfix = mpc.gen2(:,PFIX) == 1;
    idxQfix = mpc.gen2(:,QFIX) == 1;
    
    idxTrip = ~cs.activeGenerators(:,i);
    
    idxPvar = and(~idxTrip,~idxPfix);
    idxQvar = and(~idxTrip,~idxQfix);
    idxActiveGenerators = cs.activeGenerators(:,i);
    idxTrippedGenerators = ~cs.activeGenerators(:,i);
    nActiveGenerators = cs.nActiveGenerators(i);
    
    if i > 1
        Pg(idxPvar,i) = x(iPg);
        Pg(idxPfix,i) = Pg(idxPfix,1); % store fixed values in all contingencies
        Pg(idxTrip,i) = NaN;
        PgL(idxPvar,i) = Lambda.lower(iPg);
        PgU(idxPvar,i) = Lambda.upper(iPg);
        Qg(idxQvar,i) = x(iQg);
        Qg(idxQfix,i) = Qg(idxQfix,1); % store fixed values in all contingencies
        Qg(idxTrip,i) = NaN;
        QgL(idxQvar,i) = Lambda.lower(iQg);
        QgU(idxQvar,i) = Lambda.upper(iQg);
    else % base case

        Pg(:,i) = x(iPg);
        PgL(:,i) = Lambda.lower(iPg);
        PgU(:,i) = Lambda.upper(iPg);
        Qg(:,i) = x(iQg);
        QgL(:,i) = Lambda.lower(iQg);
        QgU(:,i) = Lambda.upper(iQg);
    end
        
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
Pg = [mpc.order.bus.i2e(mpc.gen(:,GEN_BUS)) Pg*mpc.baseMVA];
Qg = [mpc.order.bus.i2e(mpc.gen(:,GEN_BUS)) Qg*mpc.baseMVA];
busExt = mpc.order.bus.i2e(mpc.gen(:,GEN_BUS));
PgU = [busExt PgU];
PgL = [busExt PgL];
QgU = [busExt QgU];
QgL = [busExt QgL];


s = mpc.order.gen.e2i;
% sort according to external indexing
Pg = Pg(s,:);
Qg = Qg(s,:);
PgU = PgU(s,:);
PgL = PgL(s,:);
QgU = QgU(s,:);
QgL = QgL(s,:);

% put bus nr in Va
Va = [mpc.order.bus.i2e(mpc.bus(:,BUS_I)) Va*180/pi];
Vm = [mpc.order.bus.i2e(mpc.bus(:,BUS_I)) Vm];
busExt = mpc.order.bus.i2e(mpc.bus(:,BUS_I));

VmU = [busExt VmU];
VmL = [busExt VmL];
VaU = [busExt VaU];
VaL = [busExt VaL];

Sflow = [(1:size(branch,1))' branch(:,[F_BUS T_BUS RATE_A]) Sflow];
Sflow_lam = [(1:size(branch,1))' branch(:,[F_BUS T_BUS RATE_A]) Sflow_lam];

% set multipliers below threshold to 0
thrs = optns.lamdaTolerance;
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

varnames_pg = {'BUS'};
varnames_qg = {'BUS'};
varnames_va = {'BUS'};
varnames_vm = {'BUS'};
varnames_S = {'BRANCH','FROM','TO','RATE_A'};
for i=1:nc
	if i == 1
		stringIdx = '';
	else
		stringIdx = num2str(i);
	end
	
    varnames_pg{i+1} = ['PG' stringIdx];
    varnames_qg{i+1} = ['QG' stringIdx];
    varnames_va{i+1} = ['VA' stringIdx];
    varnames_vm{i+1} = ['VM' stringIdx];
    varnames_S{2*i+3} = ['SF' stringIdx];
    varnames_S{2*i+4} = ['ST' stringIdx];
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

S_table = array2table(Sflow,'VariableNames',varnames_S);
S_lam_table = array2table(Sflow_lam,'VariableNames',varnames_S);

tab = struct();
[tab.Pg,tab.Qg,tab.Va,tab.Vm, tab.S, tab.Slam] = deal(Pg_table,Qg_table,Va_table,Vm_table, S_table, S_lam_table);
[tab.PgUlam,tab.PgLlam,tab.QgUlam,tab.QgLlam,tab.VmUlam,tab.VmLlam,tab.VaUlam,tab.VaLlam] = ...
    deal(Pg_lamU_table,Pg_lamL_table,Qg_lamU_table,Qg_lamL_table,Vm_lamU_table,Vm_lamL_table,Va_lamU_table,Va_lamL_table);
%display(Pg_table);

%[Vm mpc.bus(:,VMIN) mpc.bus(:,VMAX) Lambda.upper(vv.i1.Vm:vv.iN.Vm)]
%runopf(mpc);

%% STORE BASE RESULTS IN MPC STRUCT

V = Vm(:,2) .* exp(1j * pi/180*Va(:,2));

% store base case voltages
bus(:,VM) = Vm(:,2);
bus(:,VA) = Va(:,2);
% store base case generation
gen(:,PG) = Pg(:,2);
gen(:,QG) = Qg(:,2);

% copy bus voltages back to gen matrix
gen(:, VG) = bus(gen(:, GEN_BUS), VM);

%% compute branch flows
Sf = V(branch(:, F_BUS)) .* conj(Yf * V);  %% cplx pwr at "from" bus, p.u.
Sflow = V(branch(:, T_BUS)) .* conj(Yt * V);  %% cplx pwr at "to" bus, p.u.
branch(:, PF) = real(Sf) * baseMVA;
branch(:, QF) = imag(Sf) * baseMVA;
branch(:, PT) = real(Sflow) * baseMVA;
branch(:, QT) = imag(Sflow) * baseMVA;

%% multipliers for branch flows
muSf = zeros(nl, 1);
muSt = zeros(nl, 1);

%% update Lagrange multipliers
bus(:, MU_VMAX)  = Lambda.upper(vv.i1.Vm:vv.iN.Vm);
bus(:, MU_VMIN)  = Lambda.lower(vv.i1.Vm:vv.iN.Vm);
gen(:, MU_PMAX)  = Lambda.upper(vv.i1.Pg:vv.iN.Pg) / baseMVA;
gen(:, MU_PMIN)  = Lambda.lower(vv.i1.Pg:vv.iN.Pg) / baseMVA;
gen(:, MU_QMAX)  = Lambda.upper(vv.i1.Qg:vv.iN.Qg) / baseMVA;
gen(:, MU_QMIN)  = Lambda.lower(vv.i1.Qg:vv.iN.Qg) / baseMVA;
bus(:, LAM_P)    = Lambda.eqnonlin(nn.i1.Pmis:nn.iN.Pmis) / baseMVA;
bus(:, LAM_Q)    = Lambda.eqnonlin(nn.i1.Qmis:nn.iN.Qmis) / baseMVA;
branch(:, MU_SF) = muSf / baseMVA;
branch(:, MU_ST) = muSt / baseMVA;

results = mpc;
[results.bus, results.branch, results.gen, results.gen2, results.bus2, ...
    results.om] = ...
        deal(bus, branch, gen, gen2, bus2, om);
% [results.bus, results.branch, results.gen, ...
%     results.om, results.x, results.mu, results.f] = ...
%         deal(bus, branch, gen, om, x, mu, f);
%branch(:, MU_SF) = muSf / baseMVA;
%branch(:, MU_ST) = muSt / baseMVA;

% %% package up results
% nlnN = getN(om, 'nln');
% nlt = length(ilt);
% ngt = length(igt);
% nbx = length(ibx);
% 
% %% extract multipliers for nonlinear constraints
% kl = find(Lambda.eqnonlin < 0);
% ku = find(Lambda.eqnonlin > 0);
% nl_mu_l = zeros(nlnN, 1);
% nl_mu_u = [zeros(2*nb, 1); muSf; muSt];
% nl_mu_l(kl) = -Lambda.eqnonlin(kl);
% nl_mu_u(ku) =  Lambda.eqnonlin(ku);
% 
% %% extract multipliers for linear constraints
% kl = find(Lambda.eqlin < 0);
% ku = find(Lambda.eqlin > 0);
% 
% mu_l = zeros(size(u));
% mu_l(ieq(kl)) = -Lambda.eqlin(kl);
% mu_l(igt) = Lambda.ineqlin(nlt+(1:ngt));
% mu_l(ibx) = Lambda.ineqlin(nlt+ngt+nbx+(1:nbx));
% 
% mu_u = zeros(size(u));
% mu_u(ieq(ku)) = Lambda.eqlin(ku);
% mu_u(ilt) = Lambda.ineqlin(1:nlt);
% mu_u(ibx) = Lambda.ineqlin(nlt+ngt+(1:nbx));
% 
% mu = struct( ...
%   'var', struct('l', Lambda.lower, 'u', Lambda.upper), ...
%   'nln', struct('l', nl_mu_l, 'u', nl_mu_u), ...
%   'lin', struct('l', mu_l, 'u', mu_u) );
% 
% results = mpc;
% [results.bus, results.branch, results.gen, ...
%     results.om, results.x, results.mu, results.f] = ...
%         deal(bus, branch, gen, om, x, mu, f);
% 
% pimul = [ ...
%   results.mu.nln.l - results.mu.nln.u;
%   results.mu.lin.l - results.mu.lin.u;
%   -ones(ny>0, 1);
%   results.mu.var.l - results.mu.var.u;
% ];
% raw = struct('xr', x, 'pimul', pimul, 'info', info, 'output', Output);
% 
% %end
% %% angle limit constraint multipliers
% if ~sdp && ll.N.ang > 0
%     iang = userdata(om, 'iang');
%     mpc.branch(iang, MU_ANGMIN) = mpc.mu.lin.l(ll.i1.ang:ll.iN.ang) * pi/180;
%     mpc.branch(iang, MU_ANGMAX) = mpc.mu.lin.u(ll.i1.ang:ll.iN.ang) * pi/180;
% end
%     if ~sdp
%         %% assign values and limit shadow prices for variables
%         om_var_order = get(om, 'var', 'order');
%         for k = 1:length(om_var_order)
%             name = om_var_order(k).name;
%             if getN(om, 'var', name)
%                 idx = vv.i1.(name):vv.iN.(name);
%                 mpc.var.val.(name) = mpc.x(idx);
%                 mpc.var.mu.l.(name) = mpc.mu.var.l(idx);
%                 mpc.var.mu.u.(name) = mpc.mu.var.u(idx);
%             end
%         end
%         
%         %% assign shadow prices for linear constraints
%         om_lin_order = get(om, 'lin', 'order');
%         for k = 1:length(om_lin_order)
%             name = om_lin_order(k).name;
%             if getN(om, 'lin', name)
%                 idx = ll.i1.(name):ll.iN.(name);
%                 mpc.lin.mu.l.(name) = mpc.mu.lin.l(idx);
%                 mpc.lin.mu.u.(name) = mpc.mu.lin.u(idx);
%             end
%         end
%         
%         %% assign shadow prices for nonlinear constraints
%         if ~dc
%             om_nln_order = get(om, 'nln', 'order');
%             for k = 1:length(om_nln_order)
%                 name = om_nln_order(k).name;
%                 if getN(om, 'nln', name)
%                     idx = nn.i1.(name):nn.iN.(name);
%                     mpc.nln.mu.l.(name) = mpc.mu.nln.l(idx);
%                     mpc.nln.mu.u.(name) = mpc.mu.nln.u(idx);
%                 end
%             end
%         end
%         
%         %% assign values for components of user cost
%         om_cost_order = get(om, 'cost', 'order');
%         for k = 1:length(om_cost_order)
%             name = om_cost_order(k).name;
%             if getN(om, 'cost', name)
%                 mpc.cost.(name) = compute_cost(om, mpc.x, name);
%             end
%         end
%         
%         %% if single-block PWL costs were converted to POLY, insert dummy y into x
%         %% Note: The "y" portion of x will be nonsense, but everything should at
%         %%       least be in the expected locations.
%         pwl1 = userdata(om, 'pwl1');
%         if ~isempty(pwl1) && ~strcmp(alg, 'TRALM') && ~(strcmp(alg, 'PDIPM') && mpopt.pdipm.step_control)
%             %% get indexing
%             vv = get_idx(om);
%             if dc
%                 nx = vv.iN.Pg;
%             else
%                 nx = vv.iN.Qg;
%             end
%             y = zeros(length(pwl1), 1);
%             raw.xr = [ raw.xr(1:nx); y; raw.xr(nx+1:end)];
%             mpc.x = [ mpc.x(1:nx); y; mpc.x(nx+1:end)];
%         end
%         
%     end