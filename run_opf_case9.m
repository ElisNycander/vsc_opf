close all;
clear;
define_constants;
vscopf_define_constants;

%% OPTIONS
optns = struct();

optns.outputFile = '';

optns.caseFile = 'case9';
eval(['mpc=' optns.caseFile ';']);
optns.contingencyFile = 'case9_contingencies';
optns.stabilityMargin = 0.1;

optns.hessian = 0;
optns.verify = 1; % verify_solutions

optns.lambdaTolerance = 1e-4; % round smaller lambda to 0 for tables



optns.hessian = 0;
optns.verify = 1; % verify_solutions
optns.verifyPQconversion = 0; % convert PV buses with generators at Q-limits to PQ in power flow

optns.lambdaTolerance = 1e-4; % round smaller lambda to 0 for tables

%optns.useInitialPF = 1; % initial power flow must be solvable

%optns.useOlaussonScenarios = 1; % Scenario from Olausson (2015)
%optns.setPenetration = 1; % Manually set penetration ratio of system for base case to this value
%optns.penetrationLevel = 0.5;
%optns.replaceGeneration = 1; % replace wind with synchronous generation in base case
%optns.generationReplacementTol = 5; % in MW
%optns.windScenario = 'C1'; % 
%slackFactor = 0.9; % scale down generation, to get positive production at slack bus in base case

optns.gen.optimizeBaseP = 1; 
optns.gen.fixBaseWind = 0; % fix curtailable P for base case (only when optimizBbaseP) - NOT IMPLEMENTED
optns.gen.usePQConstraints = 0;

optns.transferCorridors = {
    [1 4;3 6;2 8]
};
optns.externalBuses = [
];

optns.saveFigures = 1;
optns.saveData = 1;
optns.caseName = 'case9_default';

%% generator options
% extra generators           
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	           
optns.gen.extra = [
     4   50 0   50       -50       1   100     1       inf     0
     6   50 0   50       -50       1   100     1       inf     0
     8   50 0   50       -50       1   100     1       inf     0
];

optns.gen.extra(:,[QMIN QMAX]) = 0;

%% the active power of generators can either be:
% 1 Variable - can vary freely for all contingencies
% 2 Fixed - set by base case power flow
% 3 Curtailable - can be curtailed relative to base case in contingencies
% Note: Variable is default
optns.gen.fixedP = [1];
optns.gen.curtailableP = [4 5 6];
optns.gen.variableP = [];
% include base case P-values in optimization or not, may take given values as "market outcome"

%optns.gen.maxPg = [1:6]; % generators for which to max production (must be fixed)
%optns.gen.maxPgLim = [3000];

% wind power scenarios - stored in columns, with as many rows as there are
% curtailable generators, or just one row, in which case the same scenario
% is applied to all curtailable generators
optns.gen.windScenarios = [
    350
];
% probabilities for wind scenarios, empty means all scenarios equally
% likely
optns.gen.windProbabilities = [
];
% activate/deactivate given wind power scenarios, if active the scenarios
% will be constructed as [wind scenarios x contingencies]
optns.gen.useWindScenarios = 1;

optns.gen.pqFactor = 0.01*ones(6,1);
%optns.gen.pqFactor = [2; 0.05; ones(1,1); 1/10*ones(3,1)];
%optns.gen.pqFactor = [0; zeros(1,1); 0.01; zeros(3,1)];

%% branch options
optns.branch.limit = 1; % turn on/off branch limits 
optns.branch.rateA = [ % branch limits, 0 means line is unconstrained
    130*ones(3,1);
    150;
    170*ones(5,1);
];

optns.branch.duplicate = []; % duplicate these branches

%% load increase
%optns.load.loadIncreaseArea = 1; % areas where to increase load
optns.bus.loadIncrease = [5 7 9]; % buses with load increase for contingencies


%% matpower options
optns.mpopt = mpoption();

optns.mpopt.pf.enforce_q_lims = 1;
optns.mpopt.opf.flow_lim = 'S';

%% solver options
foptions = optimoptions('fmincon','Algorithm','interior-point','GradObj','on','GradConstr','on');

foptions.Display = 'off'; % off, testing, iter
foptions.TolCon = 1e-10; % high value may give non-zero lagrange multipliers also for inactive constraints
foptions.TolFun = 1e0;
foptions.TolX = 1e-10;
foptions.MaxIter = 5000;

%% CHECK OPTIONS
optns = check_opf_options(mpc,optns);

%% SETUP MATPOWER CASE
mpc = setup_mpc(mpc,optns);


%% RUN INITIAL POWER FLOW FOR BASE CASE
% do pfs quietly
optns.mpopt.out.all = 0;
optns.mpopt.verbose = 0;
mpci = runpf(mpc,optns.mpopt);
assert(mpci.success == 1,'Initial power flow was unsuccessful');
% copy base case values, for initial guess
mpc.gen(:,[QG PG]) = mpci.gen(:,[QG PG]);
mpc.bus(:,[VM VA]) = mpci.bus(:,[VM VA]);


%% convert to internal indexing
mpc = ext2int(mpc);

%% CONTINGENCIES
mpc = setup_contingencies(mpc,optns); 




%% SETUP OPTIMIZATION OBJECT
om = setup_opf(mpc,optns);

disp('hej');
% objective function
f_fcn = @(x)vscopf_f_minCurtail(x, om);
% constraint function
g_fcn = @(x)vscopf_g(x, om, optns.mpopt);
if optns.hessian % hessian function
    h_fcn = @(x,lambda)vscopf_h(x,lambda,om,optns.mpopt);
    foptions = optimoptions(foptions,'Hessian','on','HessFcn',h_fcn);
end


[x0,LB,UB] = getv(om);
[A,L,U] = linear_constraints(om);

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
  fmincon(f_fcn, x0, A, U, [], [], LB, UB, g_fcn, foptions);


et = etime(clock,t0);
fprintf(['Time: %0.3f s'],et)
fprintf(['\nObjective function: %0.1f'],f*mpc.baseMVA);
%[x0 x Lambda.lower Lambda.upper LB UB]

[results,restab] = get_opf_results(om,x,Lambda,optns);
[results.et, results.f, results.success] = deal(et,f,success);

if optns.verify == 1
    success = verify_solutions(restab,om,optns);
    %fprintf(['\nSolutions verified: %i'],success);
    fprintf(['\nSolutions verified: %i'],sum(success)==mpc.contingencies.N);
end
%fd = fopen('output.txt','at');
%printpf(results,optns.fileID);'

% vv = get_idx(om);
% tab.Qg
% 180/pi*x(vv.i1.Va1:vv.iN.Va1)
 restab.Vm
 restab.Pg
% restab.Beta
% restab.Curtail
% restab.ExpCurtail
 restab.Qg
% restab.S
% restab.Slam
%restab.lamInfo
restab.PQ

success

[h,g] = g_fcn(x);
g_dev = sum(abs(g))
h_dev = min(h)

plot_results(restab,optns);