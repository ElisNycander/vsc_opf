close all;
clear;
define_constants;
vscopf_define_constants;

%% OPTIONS
optns = struct();

optns.outputFile = '';

optns.caseFile = 'case9';
optns.contingencyFile = 'case9_contingencies';
optns.stabilityMargin = 0.2;

%% generator options
optns.gen.fixPg = [3]; % generators which is fixed for stressed cases
optns.gen.fixQg = [];

optns.gen.maxPg = [3]; % generators for which to max production (must be fixed)
optns.gen.maxPgLim = [3000];

optns.branch.limit = 0; % turn on/off branch limits 
optns.branch.rateA = []; % branch limits

optns.bus.loadIncrease = 1:9; % areas where to increase load
%% matpower options
optns.mpopt = mpoption();

optns.mpopt.pf.enforce_q_lims = 1;
optns.mpopt.opf.flow_lim = 'P';

optns.mpopt.opf.tol = 1e-15;
%% CHECK OPTIONS
optns = check_opf_options(optns);

%% SETUP MATPOWER CASE
mpc = setup_mpc(optns);

%% CONTINGENCIES
mpc = setup_contingencies(mpc,optns);

% convert to internal indexing
mpc = ext2int(mpc);

%% SETUPT OPTIMIZATION OBJECT
om = setup_opf(mpc,optns);

% constraint function
g_fcn = @(x)vscopf_g(x, om, optns.mpopt);

% hessian function
h_fcn = @(x,lambda)vscopf_h(x,lambda,om,optns.mpopt);
            
%% Find solution for all contingencies
[nc,vv,baseMVA] = deal(mpc.contingencies.N,get_idx(om),mpc.baseMVA);

nvb = sum(mpc.bus2(:,LOAD_INCREASE_AREA));
nvb_idx = mpc.bus2(:,LOAD_INCREASE_AREA) == 1;
ng = size(mpc.gen,1);
nb = size(mpc.bus,1);
nPfix = sum(mpc.gen2(:,PFIX));
Pfix_idx = find(mpc.gen2(:,PFIX));
nVar = 2*nb + 2*ng;
nVarTot = (nc+1)*nVar - nc*nPfix;
nMis = 2*(nc+1)*nb;

% Matpower jacobian:    Vm, Va       Qg      Pg base, Pg cont

mdg = sparse(nMis,nVarTot);

for i=1:nc+1
    
if i == 1
    sidx = '';
else
    sidx = num2str(i-1);
end

    % get index ranges
        iVa = vv.i1.(['Va' sidx]):vv.iN.(['Va' sidx]);
        iVm = vv.i1.(['Vm' sidx]):vv.iN.(['Vm' sidx]);
        iPg = vv.i1.(['Pg' sidx]):vv.iN.(['Pg' sidx]);
        iQg = vv.i1.(['Qg' sidx]):vv.iN.(['Qg' sidx]);
        
        iPgFix = vv.i1.Pg-1+Pfix_idx;
        
        impc = get_mpc(om);
        impc.bus(nvb_idx,[PD QD]) = impc.contingencies.load(1+(i-1)*nvb:i*nvb,:);
        
        impc = runpf(impc,optns.mpopt); % run pf
        
        % store results in x
        % note: for case9 internal indexing and external indexing is the
        % same, no need for converting
        x(iVa) = pi/180*impc.order.int.bus(:,VA);
        x(iVm) = impc.order.int.bus(:,VM);
        if i==1
            x(iPg) = impc.order.int.gen(:,PG)/baseMVA;
        else
            x(iPg) = impc.order.int.gen(setdiff(1:ng,impc.order.gen.e2i(optns.gen.fixPg)),PG)/baseMVA;
        end
        x(iQg) = impc.order.int.gen(:,QG)/baseMVA;
        
        % matpower jacobian
        [Ybus,Yf,Yt] = makeYbus(baseMVA,impc.bus,impc.branch);
        if i==1 % base case
            [~,ig,~,idg] = opf_consfcn(x',om,Ybus,Yf,Yt,optns.mpopt,[]);
        else
            % x-values for building matpower jacobian
            xtmp(vv.i1.Va:vv.iN.Va) = pi/180*impc.order.int.bus(:,VA);
            xtmp(vv.i1.Vm:vv.iN.Vm) = impc.order.int.bus(:,VM);
            xtmp(vv.i1.Pg:vv.iN.Pg) = impc.order.int.gen(:,PG)/baseMVA;
            xtmp(vv.i1.Qg:vv.iN.Qg) = impc.order.int.gen(:,QG)/baseMVA;

            [~,ig,~,idg] = opf_consfcn(xtmp',om,Ybus,Yf,Yt,optns.mpopt,[]);
        end
        idg = idg'; % transpose
        if i == 1
            mdg(1:2*nb,1:nVar) = idg;
        else
            % remove and save columns corresponding to fixed Pg
            cfix = idg(:,iPgFix);
            idg(:,iPgFix) = [];
            
            mdg(1+2*(i-1)*nb:2*i*nb, 1+(i-1)*nVar:i*nVar-(i-1)*nPfix) = idg;
            mdg(1+2*(i-1)*nb:2*i*nb, iPgFix) = cfix;
        end
end
mdg = mdg';


[~,g,~,dg] = g_fcn(x'); % evaluate constraints and jacobian


mdg = full(mdg);
dg = full(dg);

% compare jacobians
[i,j] = ind2sub(size(mdg),find(mdg ~= dg));