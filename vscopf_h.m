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

%% unpack data
lim_type = upper(mpopt.opf.flow_lim(1));
mpc = get_mpc(om);
[baseMVA, bus, bus2, gen, gen2, branch,cs] = ...
    deal(mpc.baseMVA, mpc.bus, mpc.bus2, mpc.gen, mpc.gen2, mpc.branch, mpc.contingencies);
[Ybus,Yf,Yt,load] = deal(cs.Ybus,cs.Yf, cs.Yt, cs.load);
[vv,~,nn] = get_idx(om);

%% unpack needed parameters
nb = size(bus, 1);          %% number of buses
%nl = size(branch, 1);       %% number of branches
%ng = size(gen, 1);          %% number of dispatchable injections
nxyz = length(x);           %% total number of control vars of all types
nc = cs.N; % includes contingencies


if isfield(vv.i1,'Vsp')
    useVoltageControl = 1;
else
    useVoltageControl = 0;
end

%% find constrained lines - preparation for hessian of branch flow constriants
%idxConstrainedLines = branch(:,RATE_A) ~= 0; %% logical index for constrained lines
%nConstrainedLines = sum(idxConstrainedLines);           %% number of constrained lines
idxConstrainedLines = cs.constrainedActiveLines(:,1);
nConstrainedLines = cs.nConstrainedActiveLines(1);

%% containers for hessian

%find buses in load increas area
varload_idx = find(bus2(:,LOAD_INCREASE_AREA));
busVarN = length(varload_idx);
% pfix_idx = find(gen2(:,PFIX));
% pfix_N = sum(gen2(:,PFIX));
%idxPvar = gen2(:,PFIX) == 0;
%idxQvar = gen2(:,QFIX) == 0;

% get indices
nlam = length(lambda.eqnonlin) / 2 / nc;

d2G = zeros(nxyz);
d2H = zeros(nxyz);

counter = 0; % counter for entries in d2G, d2H
ycounter = 0; % counter for rows in Yf, Yt
hcounter = 0; % counter for rows in lambda, i.e. active lines with constraints
for i=1:nc %% loop over all cases
       
    if i == 1 % base case
        sidx = '';
    else % contingency
        sidx = num2str(i);
    end
    
    iPmis = nn.i1.(['Pmis' sidx]):nn.iN.(['Pmis' sidx]);
    iQmis = nn.i1.(['Qmis' sidx]):nn.iN.(['Qmis' sidx]);
    if useVoltageControl
       iCp = nn.i1.(['Cp' sidx]):nn.iN.(['Cp' sidx]);
       iCn = nn.i1.(['Cn' sidx]):nn.iN.(['Cn' sidx]);
    %   iVbal = nn.i1.(['Vbal' sidx]):nn.iN.(['Vbal' sidx]);
       iQg = vv.i1.(['Qg' sidx]):vv.iN.(['Qg' sidx]);
       iVp = vv.i1.(['Vp' sidx]):vv.iN.(['Vp' sidx]);
       iVn = vv.i1.(['Vn' sidx]):vv.iN.(['Vn' sidx]);
    end
    
    iVa = vv.i1.(['Va' sidx]):vv.iN.(['Va' sidx]);
    iVm = vv.i1.(['Vm' sidx]):vv.iN.(['Vm' sidx]);
    %iPg = vv.i1.(['Pg' sidx]):vv.iN.(['Pg' sidx]);
    %iQg = vv.i1.(['Qg' sidx]):vv.iN.(['Qg' sidx]);
    nBranchLimits = nn.N.(['Sf' sidx]) + nn.N.(['St' sidx]);
    nVar = vv.N.(['Pg' sidx]) + vv.N.(['Qg' sidx]) + vv.N.(['Vm' sidx]) + vv.N.(['Va' sidx]);
    
    % NOTE: Admittance matrices for different contingencies stack in a
    % single column
    nActiveLines = cs.nActiveLines;
    iYbus = Ybus(1+(i-1)*nb:i*nb,:);
    iYf = Yf(ycounter+1:ycounter+nActiveLines,:);
    iYt = Yt(ycounter+1:ycounter+nActiveLines,:);
    % pick out lines with active constraints out of lines active for this
    % contingency
    idxYbranch = idxConstrainedLines(cs.activeLines(:,i)); % indices for Yf and Yt
	idxBranch = cs.constrainedActiveLines(:,i);
    nBranch = cs.nConstrainedActiveLines(i);
    iiYf = iYf(idxYbranch,:);
    iiYt = iYt(idxYbranch,:);
   

    %Pg = x(iPg);
    %Qg = x(iQg);
    Va = x(iVa);
    Vm = x(iVm);
    V = Vm .* exp(1j * Va);
    
    
    bus(varload_idx,[PD QD]) = load(1+(i-1)*busVarN:i*busVarN,:);
    
    %% make hessian for full set of variables
    nxtra = nVar - 2*nb;
    
    % NOTE: Power flow equations linear in Pg and Qg, hence no modification
    % to hessian necessary. Only number of filling zeros have to be adjusted 
    % to account for the number of variable generators for a particular 
    % contingency. Same is true for hessian of branch flow limits.
    
    %% ----- evaluate Hessian of power balance constraints -----
    % select lambdas for this set of power balance constraints
    
    %lamP = lambda.eqnonlin(1+(i-1)*nlam:i*nlam);
    %lamQ = lambda.eqnonlin(1+i*nlam:(i+1)*nlam);
    lamP = lambda.eqnonlin(iPmis);
    lamQ = lambda.eqnonlin(iQmis);
    
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
    nZerosG = nnz(d2G);
    if useVoltageControl
       % Educated guess: put Lagrange multipliers for Cp and Cn on places for
       % (Vp, Vn) and Q in d2G
       d2G(sub2ind(size(d2G),iQg,iVp)) = lambda.eqnonlin(iCp);
       d2G(sub2ind(size(d2G),iVp,iQg)) = lambda.eqnonlin(iCp);
       d2G(sub2ind(size(d2G),iQg,iVn)) = lambda.eqnonlin(iCn);
       d2G(sub2ind(size(d2G),iVn,iQg)) = lambda.eqnonlin(iCn);
        
       %assert(nnz(d2G))
    end
    
    %% evaluate Hessian of branch flow constraints 
    if nConstrainedLines
            muF = lambda.ineqnonlin(hcounter+1:hcounter+nBranch);
            muT = lambda.ineqnonlin(1+hcounter+nBranch:hcounter+2*nBranch);
            %muF = lambda.ineqnonlin(1+(i-1)*nmu:i*nmu);
            %muT = lambda.ineqnonlin(1+i*nmu:(i+1)*nmu);
           % NOTE: BETTER TO USE INDICES FROM nn
            %iSf = nn.i1.(['Sf' sidx]):nn.iN.(['Sf' sidx]);
            %iSt = nn.i1.(['St' sidx]):nn.iN.(['St' sidx]);
            
            % keep all branches until last assignemtn
        if lim_type == 'I'          %% square of current
            [dIf_dVa, dIf_dVm, dIt_dVa, dIt_dVm, If, It] = dIbr_dV(branch(idxBranch,:), iiYf, iiYt, V);
            [Hfaa, Hfav, Hfva, Hfvv] = d2AIbr_dV2(dIf_dVa, dIf_dVm, If, iiYf, V, muF);
            [Htaa, Htav, Htva, Htvv] = d2AIbr_dV2(dIt_dVa, dIt_dVm, It, iiYt, V, muT);
        else
            f = branch(idxBranch, F_BUS);    %% list of "from" buses
            t = branch(idxBranch, T_BUS);    %% list of "to" buses
            Cf = sparse(1:nBranch, f, ones(nBranch, 1), nBranch, nb);     %% connection matrix for line & from buses
            Ct = sparse(1:nBranch, t, ones(nBranch, 1), nBranch, nb);     %% connection matrix for line & to buses
            [dSf_dVa, dSf_dVm, dSt_dVa, dSt_dVm, Sf, St] = dSbr_dV(branch(idxBranch,:), iiYf, iiYt, V);
            if lim_type == '2'        %% square of real power
                [Hfaa, Hfav, Hfva, Hfvv] = d2ASbr_dV2(real(dSf_dVa), real(dSf_dVm), real(Sf), Cf, iiYf, V, muF);
                [Htaa, Htav, Htva, Htvv] = d2ASbr_dV2(real(dSt_dVa), real(dSt_dVm), real(St), Ct, iiYt, V, muT);
            elseif lim_type == 'P'    %% real power
                [Hfaa, Hfav, Hfva, Hfvv] = d2Sbr_dV2(Cf, iiYf, V, muF);
                [Htaa, Htav, Htva, Htvv] = d2Sbr_dV2(Ct, iiYt, V, muT);
                [Hfaa, Hfav, Hfva, Hfvv] = deal(real(Hfaa), real(Hfav), real(Hfva), real(Hfvv));
                [Htaa, Htav, Htva, Htvv] = deal(real(Htaa), real(Htav), real(Htva), real(Htvv));
            else                      %% square of apparent power
                [Hfaa, Hfav, Hfva, Hfvv] = d2ASbr_dV2(dSf_dVa, dSf_dVm, Sf, Cf, iiYf, V, muF);
                [Htaa, Htav, Htva, Htvv] = d2ASbr_dV2(dSt_dVa, dSt_dVm, St, Ct, iiYt, V, muT);
            end
        end
        d2H(1+counter:nVar+counter,1+counter:nVar+counter) = [
            [Hfaa Hfav; Hfva Hfvv] + [Htaa Htav; Htva Htvv] sparse(2*nb, nxtra);
            sparse(nxtra, 2*nb + nxtra)
            ];
        
    end
    counter = counter + nVar;
    hcounter = hcounter + nBranchLimits;
    ycounter = ycounter + nActiveLines;
end



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
