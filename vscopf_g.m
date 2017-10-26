function [h, g, dh, dg] = vscopf_g(x, om, mpopt,il)
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
[baseMVA, bus, gen, branch,cs,gen2,bus2] = ...
    deal(mpc.baseMVA, mpc.bus, mpc.gen, mpc.branch,mpc.contingencies, mpc.gen2, mpc.bus2);
[Ybus, Yf, Yt, load] = deal(cs.Ybus, cs.Yf, cs.Yt, cs.load);
vv = get_idx(om);

lim_type = mpopt.opf.flow_lim; %% branch flow limit type
%% problem dimensions
nb = size(bus, 1);          %% number of buses
nl = size(branch, 1);       %% number of branches
ng = size(gen, 1);          %% number of dispatchable injections
nxyz = length(x);           %% total number of control vars of all types
nc = cs.N;


%% find constrained lines
il = find(branch(:,RATE_A) ~= 0);
nl2 = length(il);           %% number of constrained lines

%% preallocate return variables
g = zeros(2*nb*(nc+1),1);
dg = sparse(length(g),nxyz);

if nl2 % prepare for line flow constraints
    h = zeros(2*nl*(nc+1),1);
    dh = sparse(length(h),nxyz);
    
    flow_max = branch(il, RATE_A) / baseMVA;
    flow_max(flow_max == 0) = Inf;
    if lim_type ~= 'P'        %% typically use square of flow
        flow_max = flow_max.^2;
    end
else
    h = [];
    dh = [];
end

%find buses in load increas area
varload_idx = find(bus2(:,LOAD_INCREASE_AREA));
busVarN = length(varload_idx);
pfix_idx = find(gen2(:,PFIX));
pfix_N = sum(gen2(:,PFIX));
pvar_idx = find(gen2(:,PFIX) == 0);

for i=1:nc+1     % loop over all cases, including base case
    
    % use admittance matrix for this contingency
    % NOTE: Admittance matrices for different contingencies stack in a
    % single column
    iYbus = Ybus(1+(i-1)*nb:i*nb,:);
    iYf = Yf(1+(i-1)*nb:i*nb,:);
    iYt = Yt(1+(i-1)*nb:i*nb,:);
    
    if i == 1 % base case
        sidx = '';
    else % contingency
        sidx = num2str(i-1);
    end
    
    Pg = x(vv.i1.(['Pg' sidx]):vv.iN.(['Pg' sidx]));
    Qg = x(vv.i1.(['Qg' sidx]):vv.iN.(['Qg' sidx]));
    Va = x(vv.i1.(['Va' sidx]):vv.iN.(['Va' sidx]));
    Vm = x(vv.i1.(['Vm' sidx]):vv.iN.(['Vm' sidx]));
    V = Vm .* exp(1j * Va);
    
    if i == 1 % base case
        gen(:, PG) = Pg * baseMVA;  %% active generation in MW
    else
        gen(pvar_idx,PG) = Pg * baseMVA;
    end
    gen(:, QG) = Qg * baseMVA;  %% reactive generation in MVAr
    
%     % increase load for contingencies
%     if i>1
%         bus(varload_idx,[PD QD]) = ...
%             baseload(varload_idx,:) * (1 + mpc.stabilityMargin);
%     end
    bus(varload_idx,[PD QD]) = load(1+(i-1)*busVarN:i*busVarN,:);
    
    Sbus = makeSbus(baseMVA, bus, gen, mpopt, Vm);  %% net injected power in p.u.
    
    % evaluate power flow equations
    mis = V .* conj(iYbus * V) - Sbus;
    
    g(1+(i-1)*2*nb:2*i*nb) = [real(mis); imag(mis)];
    
    if nl2 > 0 % then, the inequality constraints (branch flow limits)
        
        if lim_type == 'I'    %% current magnitude limit, |I|
            If = iYf * V;
            It = iYt * V;
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
    
    %%----- evaluate partials of constraints -----
    if nargout > 2
        
        %% BASE CASE JACOBIAN
        %% index ranges
        iVa = vv.i1.(['Va' sidx]):vv.iN.(['Va' sidx]);
        iVm = vv.i1.(['Vm' sidx]):vv.iN.(['Vm' sidx]);
        iPg = vv.i1.(['Pg' sidx]):vv.iN.(['Pg' sidx]);
        iQg = vv.i1.(['Qg' sidx]):vv.iN.(['Qg' sidx]);
        PgcN = vv.N.(['Pg' sidx]);
        
        %% compute partials of injected bus powers
        [dSbus_dVm, dSbus_dVa] = dSbus_dV(iYbus, V);           %% w.r.t. V
        [~, neg_dSd_dVm] = makeSbus(baseMVA, bus, gen, mpopt, Vm); % for voltage dependent loads
        dSbus_dVm = dSbus_dVm - neg_dSd_dVm;
        if i == 1
            neg_CgP = sparse(gen(:, GEN_BUS), 1:PgcN, -1, nb, PgcN);   %% Pbus w.r.t. Pg
        else
            neg_CgP = sparse(gen(pvar_idx, GEN_BUS), 1:PgcN, -1, nb, PgcN);
        end
        neg_CgQ = sparse(gen(:,GEN_BUS), 1:ng,-1,nb,ng); %% Qbus w.r.t. Qg
        
        %% construct Jacobian of equality (power flow) constraints and transpose it
        dg(1+2*nb*(i-1):2*nb*i, [iVa iVm iPg iQg]) = [
            real([dSbus_dVa dSbus_dVm]) neg_CgP sparse(nb, ng);  %% P mismatch w.r.t Va, Vm, Pg, Qg
            imag([dSbus_dVa dSbus_dVm]) sparse(nb, PgcN) neg_CgQ;  %% Q mismatch w.r.t Va, Vm, Pg, Qg
            ];
        
        if PgcN < ng
            dg(sub2ind(size(dg),gen(pfix_idx,GEN_BUS)+2*nb*(i-1),pfix_idx+vv.i1.Pg-1)) = - ones(1, ng-PgcN);
            % add derivative terms wrt base case P for fixed generators
            %dg(1+2*nb*(i-1):nb+2*nb*(i-1),pfix_idx+vv.i1.Pg-1) = -ones(nb,ng-PgcN);
        end
    end
end
% transpose jacobians
dg = dg';
dh = dh';
%display(num2str(looptime))