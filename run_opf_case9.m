close all;
clear;
define_constants;
vscopf_define_constants;

%% OPTIONS
optns = struct();

optns.outputFile = '';

optns.caseFile = 'case9';
optns.contingencyFile = 'case9_contingencies';
optns.stabilityMargin = 0.1;

optns.hessian = 0;
optns.verify = 1; % verify_solutions

optns.lamdaTolerance = 1e-8; % round smaller lambda to 0 for tables

%% generator options
% extra generators           
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	           
optns.gen.extra = [
    4   10 0   0       0       1   100     1       1e3     0
    6   10 0   0       0       1   100     1       1e3     0
    8   10 0   0       0       1   100     1       1e3     0
];

%% the active power of generators can either be:
% 1 Variable - can vary freely for all contingencies
% 2 Fixed - set by base case power flow
% 3 Curtailable - can be curtailed relative to base case in contingencies
% Note: Variable is default
optns.gen.fixedP = [4 5];
optns.gen.curtailableP = [6];
optns.gen.variableP = [1 2 3];
% include non-variable p in optimization or not, may take given values as "market outcome"
optns.gen.optimizeBaseP = 0; 

optns.gen.maxPg = [1:6]; % generators for which to max production (must be fixed)
optns.gen.maxPgLim = [3000];

optns.branch.limit = 1; % turn on/off branch limits 
optns.branch.rateA = [ % branch limits, 0 means line is unconstrained
    250*ones(4,1)
    250
    250*ones(4,1)
];

optns.branch.duplicate = []; % duplicate these branches

%optns.load.loadIncreaseArea = 1; % areas where to increase load
optns.bus.loadIncrease = [5 7 9]; % buses with load increase for contingencies
%% matpower options
optns.mpopt = mpoption();

optns.mpopt.pf.enforce_q_lims = 1;
optns.mpopt.opf.flow_lim = 'P';

%% solver options
foptions = optimoptions('fmincon','Algorithm','interior-point','GradObj','on','GradConstr','on');

foptions.Display = 'off'; % off, testing, iter
foptions.TolCon = 1e-10; % high value may give non-zero lagrange multipliers also for inactive constraints
foptions.TolFun = 1e0;
foptions.TolX = 1e-10;
foptions.MaxIter = 200;

%% CHECK OPTIONS
optns = check_opf_options(optns);

%% SETUP MATPOWER CASE
mpc = setup_mpc(optns);


%% RUN INITIAL POWER FLOW FOR BASE CASE

mpci = runpf(mpc,optns.mpopt);
assert(mpci.success == 1,'Initial power flow was unsuccessful');
% copy base case values, for initial guess
mpc.gen(:,[QG PG]) = mpci.gen(:,[QG PG]);
mpc.bus(:,[VM VA]) = mpci.bus(:,[VM VA]);
%% CONTINGENCIES
mpc = setup_contingencies(mpc,optns); 

%% convert to internal indexing
mpc = ext2int(mpc);


%% SETUP OPTIMIZATION OBJECT
om = setup_opf(mpc,optns);

% objective function
f_fcn = @(x)vscopf_f_minCurtail(x, om);
% constraint function
g_fcn = @(x)vscopf_g(x, om, optns.mpopt);
if optns.hessian % hessian function
    h_fcn = @(x,lambda)vscopf_h(x,lambda,om,optns.mpopt);
    foptions = optimoptions(foptions,'Hessian','on','HessFcn',h_fcn);
end


[x0,LB,UB] = getv(om);

%mpc = runpf(mpc,optns.mpopt);
% x0 = [pi/180*mpc.order.int.bus(:,VA); 
%      mpc.order.int.bus(:,VM);
%      mpc.order.int.gen(:,PG)/mpc.baseMVA;
%      mpc.order.int.gen(:,QG)/mpc.baseMVA;
%      pi/180*mpc.order.int.bus(:,VA); 
%      mpc.order.int.bus(:,VM);
%      mpc.order.int.gen(mpc.order.int.gen2(:,PFIX)==0,PG)/mpc.baseMVA;
%      mpc.order.int.gen(:,QG)/mpc.baseMVA ];
%  
%  x0(x0 < LB)
%  x0(x0 > UB)

t0 = clock();
[x, f, success, Output, Lambda] = ...
  fmincon(f_fcn, x0, [], [], [], [], LB, UB, g_fcn, foptions);
et = etime(clock,t0);
fprintf(['Time: %0.3f s'],et)
fprintf(['\nObjective function: %0.1f'],-f*mpc.baseMVA);
%[x0 x Lambda.lower Lambda.upper LB UB]

[results,tab] = get_opf_results(om,x,Lambda,optns);
[results.et, results.f, results.success] = deal(et,f,success);

if optns.verify == 1
    success = verify_solutions(tab,om,optns);
    fprintf(['\nSolutions verified: %i'],success);
end
%fd = fopen('output.txt','at');
%printpf(results,optns.fileID);'

% vv = get_idx(om);
% tab.Qg
% 180/pi*x(vv.i1.Va1:vv.iN.Va1)
tab.Vm
tab.Pg
tab.Qg
tab.S
tab.Slam

[h,g] = g_fcn(x);
g_dev = sum(abs(g))
h_dev = min(abs(h))