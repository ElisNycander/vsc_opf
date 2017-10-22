function [h, g, dh, dg] = vscopf_g(x, om ,Ybus,Yf,Yt, mpopt,il)
%function [h,g,dh,dg] = pwind_consfcn(x, om, Ybus, Yf, Yt, mpopt, il, varargin)
%OPF_CONSFCN  Evaluates nonlinear constraints and their Jacobian for OPF.
%   [H, G, DH, DG] = OPF_CONSFCN(X, OM, YBUS, YF, YT, MPOPT, IL)
%
%   Constraint evaluation function for AC optimal power flow, suitable
%   for use with MIPS or FMINCON. Computes constraint vectors and their
%   gradients.
%
%   Inputs:
%     X : optimization vector
%     OM : OPF model object
%     YBUS : bus admittance matrix
%     YF : admittance matrix for "from" end of constrained branches
%     YT : admittance matrix for "to" end of constrained branches
%     MPOPT : MATPOWER options struct
%     IL : (optional) vector of branch indices corresponding to
%          branches with flow limits (all others are assumed to be
%          unconstrained). The default is [1:nl] (all branches).
%          YF and YT contain only the rows corresponding to IL.
%
%   Outputs:
%     H  : vector of inequality constraint values (flow limits)
%          where the flow can be apparent power, real power, or
%          current, depending on the value of opf.flow_lim in MPOPT
%          (only for constrained lines), normally expressed as
%          (limit^2 - flow^2), except when opf.flow_lim == 'P',
%          in which case it is simply (limit - flow).
%     G  : vector of equality constraint values (power balances)
%     DH : (optional) inequality constraint gradients, column j is
%          gradient of H(j)
%     DG : (optional) equality constraint gradients
%
%   Examples:
%       [h, g] = opf_consfcn(x, om, Ybus, Yf, Yt, mpopt);
%       [h, g, dh, dg] = opf_consfcn(x, om, Ybus, Yf, Yt, mpopt);
%       [h, g, dh, dg] = opf_consfcn(x, om, Ybus, Yf, Yt, mpopt, il);
%
%   See also OPF_COSTFCN, OPF_HESSFCN.

tic;
%% define named indices into data matrices
define_constants;

%% unpack data
mpc = get_mpc(om);
[baseMVA, bus, gen, branch,contingencies] = deal(mpc.baseMVA, mpc.bus, mpc.gen, mpc.branch,mpc.contingencies);
vv = get_idx(om);

%% problem dimensions
nb = size(bus, 1);          %% number of buses
nl = size(branch, 1);       %% number of branches
ng = size(gen, 1);          %% number of dispatchable injections
nxyz = length(x);           %% total number of control vars of all types
nc = size(contingencies,1);

%% find constrained lines
il = find(branch(:,RATE_A) ~= 0); 
nl2 = length(il);           %% number of constrained lines

%% return variable containers
g = zeros(2*nb*(nc+1),1);

%% BASE CASE MISMATCH
%% grab Pg & Qg
Pg = x(vv.i1.Pg:vv.iN.Pg);  %% active generation in p.u.
Qg = x(vv.i1.Qg:vv.iN.Qg);  %% reactive generation in p.u.

%% put Pg & Qg back in gen
gen(:, PG) = Pg * baseMVA;  %% active generation in MW
gen(:, QG) = Qg * baseMVA;  %% reactive generation in MVAr

%% ----- evaluate constraints -----
%% reconstruct V
Va = x(vv.i1.Va:vv.iN.Va);
Vm = x(vv.i1.Vm:vv.iN.Vm);
V = Vm .* exp(1j * Va);

%% rebuild Sbus
Sbus = makeSbus(baseMVA, bus, gen, mpopt, Vm);  %% net injected power in p.u.

%% evaluate power flow equations
mis = V .* conj(Ybus * V) - Sbus;

g(1:2*nb) = [real(mis); imag(mis)];


% find buses in load increas area
bcidx = find(mpc.bus2(:,LOAD_INCREASE_AREA));
gcidx = find(mpc.gen2(:,PFIX) == 0);

Vac = zeros(nb,nc);
Vmc = zeros(size(Vac));
Vc = zeros(size(Vac));
% initialize array
if nc
    mpc_array(nc) = mpc;
end
nz = nnz(g);
for i=1:nc
    mpcc = mpc; % get base case
    si = num2str(i);
    % get contingency variables
    Pgc = x(eval(['vv.i1.Pg' si ':vv.iN.Pg' si]));
    Qgc = x(eval(['vv.i1.Qg' si ':vv.iN.Qg' si]));
    Vac(:,i) = x(eval(['vv.i1.Va' si ':vv.iN.Va' si]));
    Vmc(:,i) = x(eval(['vv.i1.Vm' si ':vv.iN.Vm' si]));
    Vc(:,i) = Vmc(:,i) .* exp(1j * Vac(:,i));
    
    % construct generation
    mpcc.gen(~gcidx,PG) = gen(~gcidx,PG); % insert values from base case
    mpcc.gen(gcidx,PG) = Pgc*baseMVA; % contingency values
    mpcc.gen(:,QG) = Qgc*baseMVA;
    % construct load
    mpcc.bus(bcidx,[PD QD]) = ...
        mpcc.bus(bcidx,[PD QD]) * (1 + mpc.stabilityMargin);
    
    Sbus_c = makeSbus(baseMVA, mpcc.bus, mpcc.gen, mpopt, Vmc(:,i));
    % note: Ybus should be reconstructed with applied contingencies
    mis_c = Vc(:,i) .* conj(Ybus * V) - Sbus_c;
    
    nzr = sum(real(mis_c) ~= 0);
    nzi = sum(imag(mis_c) ~= 0);
    
    g(2*nb+1+(i-1)*2*nb:2*nb+2*i*nb) = [real(mis_c); imag(mis_c)];
    
    mpc_array(i) = mpcc;
    assert(nnz(g) == nz + nzr + nzi,'Overwrite existing entry in constraint g');
    
end


lim_type = mpopt.opf.flow_lim;
%% then, the inequality constraints (branch flow limits)
if nl2 > 0
  h = zeros(2*nl*(nc+1),1);
  
  flow_max = branch(il, RATE_A) / baseMVA;
  flow_max(flow_max == 0) = Inf;
  if lim_type ~= 'P'        %% typically use square of flow
    flow_max = flow_max.^2;
  end
  
  % loop over all cases, including base case
  for i=1:nc+1
  
  if lim_type == 'I'    %% current magnitude limit, |I|
    If = Yf * V;
    It = Yt * V;
    h = [ If .* conj(If) - flow_max;    %% branch current limits (from bus)
          It .* conj(It) - flow_max ];  %% branch current limits (to bus)
  else
    %% compute branch power flows
    Sf = V(branch(il, F_BUS)) .* conj(Yf * V);  %% complex power injected at "from" bus (p.u.)
    St = V(branch(il, T_BUS)) .* conj(Yt * V);  %% complex power injected at "to" bus (p.u.)
    if lim_type == '2'                          %% active power limit, P squared (Pan Wei)
      h = [ real(Sf).^2 - flow_max;             %% branch real power limits (from bus)
            real(St).^2 - flow_max ];           %% branch real power limits (to bus)
    elseif lim_type == 'P'                      %% active power limit, P
      h = [ real(Sf) - flow_max;                %% branch real power limits (from bus)
            real(St) - flow_max ];              %% branch real power limits (to bus
    else                                        %% apparent power limit, |S|
      h = [ Sf .* conj(Sf) - flow_max;          %% branch apparent power limits (from bus)
            St .* conj(St) - flow_max ];        %% branch apparent power limits (to bus)
    end
  end
  end
else
  h = [];
end

%%----- evaluate partials of constraints -----
if nargout > 2
    dg = sparse(2*nb*(nc+1),nxyz);
    
    %% BASE CASE JACOBIAN
    %% index ranges
    iVa = vv.i1.Va:vv.iN.Va;
    iVm = vv.i1.Vm:vv.iN.Vm;
    iPg = vv.i1.Pg:vv.iN.Pg;
    iQg = vv.i1.Qg:vv.iN.Qg;
    
    %% compute partials of injected bus powers
    [dSbus_dVm, dSbus_dVa] = dSbus_dV(Ybus, V);           %% w.r.t. V
    [~, neg_dSd_dVm] = makeSbus(baseMVA, bus, gen, mpopt, Vm); % for voltage dependent loads
    dSbus_dVm = dSbus_dVm - neg_dSd_dVm;
    neg_Cg = sparse(gen(:, GEN_BUS), 1:ng, -1, nb, ng);   %% Pbus w.r.t. Pg
    %% Qbus w.r.t. Qg
    
    %% construct Jacobian of equality (power flow) constraints and transpose it
    %dg = sparse(2*nb, nxyz);
    dg(1:2*nb, [iVa iVm iPg iQg]) = [
        real([dSbus_dVa dSbus_dVm]) neg_Cg sparse(nb, ng);  %% P mismatch w.r.t Va, Vm, Pg, Qg
        imag([dSbus_dVa dSbus_dVm]) sparse(nb, ng) neg_Cg;  %% Q mismatch w.r.t Va, Vm, Pg, Qg
        ];
    
    
    %  dg = dg';
    
    %% CONSTRAINT JACOBIANS
    nz = nnz(dg);
    for i=1:nc
        
        % get var indices
        si = num2str(i);
        iPgc = eval(['vv.i1.Pg' si ':vv.iN.Pg' si]);
        iQgc = eval(['vv.i1.Qg' si ':vv.iN.Qg' si]);
        iVac = eval(['vv.i1.Va' si ':vv.iN.Va' si]);
        iVmc = eval(['vv.i1.Vm' si ':vv.iN.Vm' si]);
        PgcN = eval(['vv.N.Pg' si]);
        QgcN = eval(['vv.N.Qg' si]);
        
        
        %% compute partials of injected bus powers
        %% NOTE: should update Ybus if needed
        mpc = mpc_array(i);
        [dSbus_dVm, dSbus_dVa] = dSbus_dV(Ybus, Vc(:,i));           %% w.r.t. V
        [~, neg_dSd_dVm] = makeSbus(baseMVA, mpc.bus, ...
            mpc.gen, ...
            mpopt, Vm(:,i)); % for voltage dependent loads
        dSbus_dVm = dSbus_dVm - neg_dSd_dVm;
        
        neg_CgP = sparse(mpc.gen(gcidx, GEN_BUS), 1:PgcN, -1, nb, PgcN);   %% Pbus w.r.t. Pg
        neg_CgQ = sparse(mpc.gen(:,GEN_BUS), 1:ng,-1,nb,ng); %% Qbus w.r.t. Qg
        
        
        dgc = [
            real([dSbus_dVa dSbus_dVm]) neg_CgP sparse(nb, ng);  %% P mismatch w.r.t Va, Vm, Pg, Qg
            imag([dSbus_dVa dSbus_dVm]) sparse(nb, PgcN) neg_CgQ;  %% Q mismatch w.r.t Va, Vm, Pg, Qg
            ];
        nz_add = nnz(dgc);
        dg([i*2*nb+1:(i+1)*2*nb], [iVac iVmc iPgc iQgc]) = dgc;
        
        % add derivative terms wrt base case P for fixed generators
        dg(i*2*nb+1:i*2*nb+nb,find(~gcidx)+vv.i1.Pg-1) = ones(nb,sum(~gcidx));
        
        assert(nz+nz_add+nb*sum(~gcidx) == nnz(dg),'Overwriting existing term in ''dg'', indexing wrong')
        
        
    end
    dg = dg';
end

dh = [];
looptime = toc;
%display(num2str(looptime))