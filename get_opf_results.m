function [results,tab] = get_opf_results(om,x,Lambda)

%% for script
% clear;
% close all;
% load('vscopf_post_processing_data.mat');

%%
define_constants;
[vv, ll, nn] = get_idx(om);
mpc = get_mpc(om);

[contingencies,gen,bus,bus2,gen2,branch,baseMVA] = deal(mpc.contingencies,mpc.gen,mpc.bus,mpc.bus2,mpc.gen2,mpc.branch,mpc.baseMVA);

nc = contingencies.N;
ng = size(gen,1);
nb = size(bus,1);
nl = size(branch,1);

[~,Yf,Yt] = makeYbus(baseMVA,bus,branch);

%% PRINT GENERATION

% collect injections in matrix:
%     P0  P1  ... Pnc (contingencies)
% G1
% G2
% ...
% Gng 
% (generators)

Pg = zeros(ng,nc+1);
Qg = zeros(ng,nc+1);
PgU = zeros(ng,nc+1);
PgL = zeros(ng,nc+1);
QgU = zeros(ng,nc+1);
QgL = zeros(ng,nc+1);

Pg(:,1) = x(vv.i1.Pg:vv.iN.Pg);
Qg(:,1) = x(vv.i1.Qg:vv.iN.Qg);
PgL(:,1) = Lambda.lower(vv.i1.Pg:vv.iN.Pg);
PgU(:,1) = Lambda.upper(vv.i1.Pg:vv.iN.Pg);
QgL(:,1) = Lambda.lower(vv.i1.Qg:vv.iN.Qg);
QgU(:,1) = Lambda.upper(vv.i1.Qg:vv.iN.Qg);

for i=1:nc
    % indices
    pidxs = eval(['vv.i1.Pg' num2str(i) ':vv.iN.Pg' num2str(i)]);
    qidxs = eval(['vv.i1.Qg' num2str(i) ':vv.iN.Qg' num2str(i)]);
    pfix = mpc.gen2(:,PFIX) == 1;
    qfix = mpc.gen2(:,QFIX) == 1;
    
    Pg(~pfix,i+1) = x(pidxs);
    Pg(pfix,i+1) = NaN;    
    PgL(~pfix,i+1) = Lambda.lower(pidxs);
    PgU(~pfix,i+1) = Lambda.upper(pidxs);
    Qg(~qfix,i+1) = x(qidxs);
    Qg(qfix,i+1) = NaN;
    QgL(~qfix,i+1) = Lambda.lower(qidxs);
    QgU(~qfix,i+1) = Lambda.upper(qidxs);
end

% put gen bus in matrix, convert to nominal units
Pg = [mpc.order.bus.i2e(mpc.gen(:,GEN_BUS)) Pg*mpc.baseMVA];
Qg = [mpc.order.bus.i2e(mpc.gen(:,GEN_BUS)) Qg*mpc.baseMVA];

s = mpc.order.gen.e2i;
% sort according to external indexing
Pg = Pg(s,:);
Qg = Qg(s,:);
PgU = PgU(s,:);
PgL = PgL(s,:);
QgU = QgU(s,:);
QgL = QgL(s,:);

%% PRINT VOLTAGES

Vm = zeros(nb,nc+1);
Va = zeros(size(Vm));
VmU = zeros(size(Vm));
VmL = zeros(size(Vm));
VaU = zeros(size(Vm));
VaL = zeros(size(Vm));

Vm(:,1) = x(vv.i1.Vm:vv.iN.Vm);
Va(:,1) = x(vv.i1.Va:vv.iN.Va);
VmL(:,1) = Lambda.lower(vv.i1.Vm:vv.iN.Vm);
VmU(:,1) = Lambda.upper(vv.i1.Vm:vv.iN.Vm);
VaL(:,1) = Lambda.lower(vv.i1.Va:vv.iN.Va);
VaU(:,1) = Lambda.upper(vv.i1.Va:vv.iN.Va);

for i=1:nc
    vm_idx = eval(['vv.i1.Vm' num2str(i) ':vv.iN.Vm' num2str(i)]);
    va_idx = eval(['vv.i1.Va' num2str(i) ':vv.iN.Va' num2str(i)]);
    Vm(:,i+1) = x(vm_idx);
    Va(:,i+1) = x(va_idx);
    
    VmU(:,i+1) = Lambda.upper(vm_idx);
    VmL(:,i+1) = Lambda.lower(vm_idx);
    VaU(:,i+1) = Lambda.upper(va_idx);
    VaL(:,i+1) = Lambda.lower(va_idx);
    
end

% put bus nr in Va
Va = [mpc.order.bus.i2e(mpc.bus(:,BUS_I)) Va*180/pi];
Vm = [mpc.order.bus.i2e(mpc.bus(:,BUS_I)) Vm];

% sort according to external indexing
%s = mpc.order.bus.
%Va = Va(

%% make tables

varnames_pg = {'BUS','PG0'};
varnames_qg = {'BUS','QG0'};
varnames_va = {'BUS','VA0'};
varnames_vm = {'BUS','VM0'};
for i=1:nc
    varnames_pg{i+2} = ['PG' num2str(i)];
    varnames_qg{i+2} = ['QG' num2str(i)];
    varnames_va{i+2} = ['VA' num2str(i)];
    varnames_vm{i+2} = ['VM' num2str(i)];
end

Pg_table = array2table(Pg,...
    'VariableNames',varnames_pg);
Qg_table = array2table(Qg,...
    'VariableNames',varnames_qg);
Vm_table = array2table(Vm,'VariableNames',varnames_vm);
Va_table = array2table(Va,'VariableNames',varnames_va);

tab = struct();
[tab.Pg,tab.Qg,tab.Va,tab.Vm] = deal(Pg_table,Qg_table,Va_table,Vm_table);

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
St = V(branch(:, T_BUS)) .* conj(Yt * V);  %% cplx pwr at "to" bus, p.u.
branch(:, PF) = real(Sf) * baseMVA;
branch(:, QF) = imag(Sf) * baseMVA;
branch(:, PT) = real(St) * baseMVA;
branch(:, QT) = imag(St) * baseMVA;

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