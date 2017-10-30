function Lxx = vscopf_h(x, lambda, om, mpopt)

% clear;
% load('hdata.mat');
% [lambda,mpopt] = deal(Lambda,optns.mpopt);

%OPF_HESSFCN  Evaluates Hessian of Lagrangian for AC OPF.
%   LXX = OPF_HESSFCN(X, LAMBDA, COST_MULT, OM, YBUS, YF, YT, MPOPT, IL)
%
%   Hessian evaluation function for AC optimal power flow, suitable
%   for use with MIPS or FMINCON's interior-point algorithm.
%
%   Inputs:
%     X : optimization vector
%     LAMBDA (struct)
%       .eqnonlin : Lagrange multipliers on power balance equations
%       .ineqnonlin : Kuhn-Tucker multipliers on constrained branch flows
%     COST_MULT : (optional) Scale factor to be applied to the cost
%          (default = 1).
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


%%----- initialize -----
define_constants;

%% default args
cost_mult = 1;

%% unpack data
lim_type = upper(mpopt.opf.flow_lim(1));
mpc = get_mpc(om);
[baseMVA, bus, bus2, gen, gen2, branch, gencost,c] = ...
    deal(mpc.baseMVA, mpc.bus, mpc.bus2, mpc.gen, mpc.gen2, mpc.branch, mpc.gencost, mpc.contingencies);
[Ybus,Yf,Yt,load] = deal(c.Ybus,c.Yf, c.Yt, c.load);
vv = get_idx(om);

%% unpack needed parameters
nb = size(bus, 1);          %% number of buses
nl = size(branch, 1);       %% number of branches
ng = size(gen, 1);          %% number of dispatchable injections
nxyz = length(x);           %% total number of control vars of all types
nc = c.N;

%% find constrained lines
il = find(branch(:,RATE_A) ~= 0);
nl2 = length(il);           %% number of constrained lines


%% containers for hessian

%find buses in load increas area
varload_idx = find(bus2(:,LOAD_INCREASE_AREA));
busVarN = length(varload_idx);
pfix_idx = find(gen2(:,PFIX));
pfix_N = sum(gen2(:,PFIX));
pvar_idx = find(gen2(:,PFIX) == 0);

d2G = zeros(nxyz);
counter = 0;
for i=1:nc+1 %% loop over all cases
    
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
    iVa = vv.i1.(['Va' sidx]):vv.iN.(['Va' sidx]);
    iVm = vv.i1.(['Vm' sidx]):vv.iN.(['Vm' sidx]);
    iPg = vv.i1.(['Pg' sidx]):vv.iN.(['Pg' sidx]);
    iQg = vv.i1.(['Qg' sidx]):vv.iN.(['Qg' sidx]);
    PgcN = vv.N.(['Pg' sidx]);
    nVar = vv.N.(['Pg' sidx]) + vv.N.(['Qg' sidx]) + vv.N.(['Vm' sidx]) + vv.N.(['Va' sidx]);
    
    Pg = x(iPg);
    Qg = x(iQg);
    Va = x(iVa);
    Vm = x(iVm);
    
    V = Vm .* exp(1j * Va);
    
    
    if i == 1 % base case
        gen(:, PG) = Pg * baseMVA;  %% active generation in MW
    else
        gen(pvar_idx,PG) = Pg * baseMVA;
    end
    gen(:, QG) = Qg * baseMVA;  %% reactive generation in MVAr
    
    bus(varload_idx,[PD QD]) = load(1+(i-1)*busVarN:i*busVarN,:);
    
    %% select bus for present contingency
    iYbus = Ybus(1+(i-1)*nb:i*nb,:);
    
    %% make hessian for full set of variables
    nxtra = nVar - 2*nb;
    
    % NOTE: Power flow equations linear in Pg and Qg, hence no modification
    % to hessian necessary. Only number of filling zeros have to be adjusted 
    % to account for the number of variable generators for a particularc 
    % contingency
    
    %% ----- evaluate Hessian of power balance constraints -----
    nlam = length(lambda.eqnonlin) / 2 / (1+nc); % select lambdas for this set of power balance constraints
    lamP = lambda.eqnonlin(1:nlam);
    lamQ = lambda.eqnonlin((1:nlam)+nlam);
    [Gpaa, Gpav, Gpva, Gpvv] = d2Sbus_dV2(iYbus, V, lamP);
    [Gqaa, Gqav, Gqva, Gqvv] = d2Sbus_dV2(iYbus, V, lamQ);
    %% constant impedance part of ZIP loads
    diaglam = sparse(1:nb, 1:nb, lamP, nb, nb);
    Sd = makeSdzip(baseMVA, bus, mpopt);
    diagSdz = sparse(1:nb, 1:nb, Sd.z, nb, nb);
    Gpvv = Gpvv + 2 * diaglam * diagSdz;
    d2G(1+counter:nVar+counter,1+counter:nVar+counter) = [
        real([Gpaa Gpav; Gpva Gpvv]) + imag([Gqaa Gqav; Gqva Gqvv]) sparse(2*nb, nxtra);
        sparse(nxtra, 2*nb + nxtra)
        ];
    counter = counter + nVar;
    
end


%% ----- evaluate d2f -----


% %% ----- evaluate Hessian of power balance constraints -----
% nlam = length(lambda.eqnonlin) / 2;
% lamP = lambda.eqnonlin(1:nlam);
% lamQ = lambda.eqnonlin((1:nlam)+nlam);
% [Gpaa, Gpav, Gpva, Gpvv] = d2Sbus_dV2(Ybus, V, lamP);
% [Gqaa, Gqav, Gqva, Gqvv] = d2Sbus_dV2(Ybus, V, lamQ);
% %% constant impedance part of ZIP loads
% diaglam = sparse(1:nb, 1:nb, lamP, nb, nb);
% Sd = makeSdzip(baseMVA, bus, mpopt);
% diagSdz = sparse(1:nb, 1:nb, Sd.z, nb, nb);
% Gpvv = Gpvv + 2 * diaglam * diagSdz;
% d2G = [
%     real([Gpaa Gpav; Gpva Gpvv]) + imag([Gqaa Gqav; Gqva Gqvv]) sparse(2*nb, nxtra);
%     sparse(nxtra, 2*nb + nxtra)
% ];

%%----- evaluate Hessian of flow constraints -----
% nmu = length(lambda.ineqnonlin) / 2;
% if nmu
%     muF = lambda.ineqnonlin(1:nmu);
%     muT = lambda.ineqnonlin((1:nmu)+nmu);
% else    %% keep dimensions of empty matrices/vectors compatible
%     muF = zeros(0,1);   %% (required to avoid problems when using Knitro
%     muT = zeros(0,1);   %%  on cases with all lines unconstrained)
% end
% if lim_type == 'I'          %% square of current
%     [dIf_dVa, dIf_dVm, dIt_dVa, dIt_dVm, If, It] = dIbr_dV(branch(il,:), Yf, Yt, V);
%     [Hfaa, Hfav, Hfva, Hfvv] = d2AIbr_dV2(dIf_dVa, dIf_dVm, If, Yf, V, muF);
%     [Htaa, Htav, Htva, Htvv] = d2AIbr_dV2(dIt_dVa, dIt_dVm, It, Yt, V, muT);
% else
%   f = branch(il, F_BUS);    %% list of "from" buses
%   t = branch(il, T_BUS);    %% list of "to" buses
%   Cf = sparse(1:nl2, f, ones(nl2, 1), nl2, nb);     %% connection matrix for line & from buses
%   Ct = sparse(1:nl2, t, ones(nl2, 1), nl2, nb);     %% connection matrix for line & to buses
%   [dSf_dVa, dSf_dVm, dSt_dVa, dSt_dVm, Sf, St] = dSbr_dV(branch(il,:), Yf, Yt, V);
%   if lim_type == '2'        %% square of real power
%     [Hfaa, Hfav, Hfva, Hfvv] = d2ASbr_dV2(real(dSf_dVa), real(dSf_dVm), real(Sf), Cf, Yf, V, muF);
%     [Htaa, Htav, Htva, Htvv] = d2ASbr_dV2(real(dSt_dVa), real(dSt_dVm), real(St), Ct, Yt, V, muT);
%   elseif lim_type == 'P'    %% real power                                 
%     [Hfaa, Hfav, Hfva, Hfvv] = d2Sbr_dV2(Cf, Yf, V, muF);
%     [Htaa, Htav, Htva, Htvv] = d2Sbr_dV2(Ct, Yt, V, muT);
%     [Hfaa, Hfav, Hfva, Hfvv] = deal(real(Hfaa), real(Hfav), real(Hfva), real(Hfvv));
%     [Htaa, Htav, Htva, Htvv] = deal(real(Htaa), real(Htav), real(Htva), real(Htvv));
%   else                      %% square of apparent power
%     [Hfaa, Hfav, Hfva, Hfvv] = d2ASbr_dV2(dSf_dVa, dSf_dVm, Sf, Cf, Yf, V, muF);
%     [Htaa, Htav, Htva, Htvv] = d2ASbr_dV2(dSt_dVa, dSt_dVm, St, Ct, Yt, V, muT);
%   end
% end
% d2H = [
%     [Hfaa Hfav; Hfva Hfvv] + [Htaa Htav; Htva Htvv] sparse(2*nb, nxtra);
%     sparse(nxtra, 2*nb + nxtra)
% ];


%%-----  do numerical check using (central) finite differences  -----
if 0
    nx = length(x);
    step = 1e-5;
    num_d2f = sparse(nx, nx);
    num_d2G = sparse(nx, nx);
    num_d2H = sparse(nx, nx);
    for i = 1:nx
        xp = x;
        xm = x;
        xp(i) = x(i) + step/2;
        xm(i) = x(i) - step/2;
        % evaluate cost & gradients
        [fp, dfp] = opf_costfcn(xp, om);
        [fm, dfm] = opf_costfcn(xm, om);
        % evaluate constraints & gradients
        [Hp, Gp, dHp, dGp] = opf_consfcn(xp, om, Ybus, Yf, Yt, mpopt, il);
        [Hm, Gm, dHm, dGm] = opf_consfcn(xm, om, Ybus, Yf, Yt, mpopt, il);
        num_d2f(:, i) = cost_mult * (dfp - dfm) / step;
        num_d2G(:, i) = (dGp - dGm) * lambda.eqnonlin   / step;
        num_d2H(:, i) = (dHp - dHm) * lambda.ineqnonlin / step;
    end
    d2f_err = full(max(max(abs(d2f - num_d2f))));
    d2G_err = full(max(max(abs(d2G - num_d2G))));
    d2H_err = full(max(max(abs(d2H - num_d2H))));
    if d2f_err > 1e-6
        fprintf('Max difference in d2f: %g\n', d2f_err);
    end
    if d2G_err > 1e-5
        fprintf('Max difference in d2G: %g\n', d2G_err);
    end
    if d2H_err > 1e-6
        fprintf('Max difference in d2H: %g\n', d2H_err);
    end
end
d2H = zeros(nxyz);
d2f = zeros(nxyz);

Lxx = d2f + d2G + d2H;
