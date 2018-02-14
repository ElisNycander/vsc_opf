close all;
clear;
define_constants;
vscopf_define_constants;

%% OPTIONS
optns = struct();

optns.outputFile = '';

optns.caseFile = 'Nordic32';
eval(['mpc=' optns.caseFile ';']); % get case

optns.contingencyFile = 'n32_contingencies';
%optns.stabilityMargin = 0.1;

optns.hessian = 0;
optns.verify = 1; % verify_solutions

optns.lambdaTolerance = 1e-4; % round smaller lambda to 0 for tables

%% generator options
% extra generators           
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	           
optns.gen.extra = [
    
%% At nodes where there is generation currently
%     4021   100 0   0       0       1   100     1       1e3     0
%     4012   100 0   0       0       1   100     1       1e3     0
%     4051   100 0   0       0       1   100     1       1e3     0
%     4047   100 0   0       0       1   100     1       1e3     0   
    
%% At nodes without generation
    1011   100 0   0       0       1   100     1       1e3     0
    2031   100 0   0       0       1   100     1       1e3     0
    41     100 0   0       0       1   100     1       1e3     0
    1041   100 0   0       0       1   100     1       1e3     0  
    
%    6   50 0   0       0       1   100     1       1e3     0
%    8   50 0   0       0       1   100     1       1e3     0
];

% for PQ capability
optns.gen.pqFactor = ones(size(mpc.gen,1)+size(optns.gen.extra,1),1);

    
%mpc.gen(11,PG) = mpc.gen(11,PG)-50;
%mpc.gen(10,PG) = mpc.gen(10,PG)-50;
%optns.gen.compBuses = [4021 4012]; % generators at these buses will have their PG decreased to compensate for extra generation in base case

%% the active power of generators can either be:
% 1 Variable - can vary freely for all contingencies
% 2 Fixed - set by base case power flow
% 3 Curtailable - can be curtailed relative to base case in contingencies
% Note: Variable is default
optns.gen.fixedP = [15 16 17 18];
optns.gen.curtailableP = [24 25 26 27];
optns.gen.variableP = [];
% include non-variable p in optimization or not, may take given values as "market outcome"
optns.gen.optimizeBaseP = 0; 


% wind power scenarios - stored in columns, with as many rows as there are
% curtailable generators, or just one row, in which case the same scenario
% is applied to all curtailable generators
optns.gen.windScenarios = [
    2000
];
% probabilities for wind scenarios, empty means all scenarios equally
% likely
optns.gen.windProbabilities = [
];
% activate/deactivate given wind power scenarios, if active the scenarios
% will be constructed as [wind scenarios x contingencies]
optns.gen.useWindScenarios = 1;

%% branch options
optns.branch.limit = 1; % turn on/off branch limits 
optns.branch.rateA = [ % branch limits, 0 means line is unconstrained, empty means default limit are used
];

% determine line rating from line voltage, see  
% Van Cutsem (2013) Description, Modeling, and Simulation Results of a Test
% System for Voltage Stability Analysis
rateFrom = mpc.branch(:,1)-rem(mpc.branch(:,1),1000);
rateTo = mpc.branch(:,2)-rem(mpc.branch(:,2),1000);
rateA = zeros(size(rateFrom));

rateA(rateFrom == rateTo) = rateFrom(rateFrom == rateTo);
% transmission lines
rateA(rateA == 1000) = 350;
rateA(rateA == 2000) = 500;
rateA(rateA == 4000) = 1400;
% step up transformers
rateA( and( rateA == 0,mpc.branch(:,RATE_A)==0 ) ) = [ ...
   1250
   1250
   833.3
   1000
   1000
   1000
   1000
   833.3
];
rateA( mpc.branch(:,RATE_A) ~= 0) = mpc.branch( mpc.branch(:,RATE_A) ~= 0,RATE_A );
%optns.branch.rateA = rateA;


optns.branch.duplicate = []; % duplicate these branches

%% load increase
%optns.load.loadIncreaseArea = 1; % areas where to increase load
optns.bus.loadIncrease = []; % buses with load increase for contingencies


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
foptions.MaxIter = 1000;

%% CHECK OPTIONS
optns = check_opf_options(optns);

%% SETUP MATPOWER CASE
mpc = setup_mpc(mpc,optns);

% set min and max values for PG (not present in N32)
mpc.gen(:,[PMAX PMIN]) = [mpc.gen(:,MBASE) zeros(size(mpc.gen,1),1)];

%% RUN INITIAL POWER FLOW FOR BASE CASE
% do pfs quietly
optns.mpopt.out.all = 0;
optns.mpopt.verbose = 0;
mpci = runpf(mpc,optns.mpopt);
assert(mpci.success == 1,'Initial power flow was unsuccessful');
% copy base case values, for initial guess
mpc.gen(:,[QG PG]) = mpci.gen(:,[QG PG]); % note: i2e
mpc.bus(:,[VM VA]) = mpci.bus(:,[VM VA]);

save('mpc_basecase.mat','mpci');

%% convert to internal indexing
mpc = ext2int(mpc);

%% CONTINGENCIES
mpc = setup_contingencies(mpc,optns); 



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
fprintf(['\nObjective function: %0.1f'],f*mpc.baseMVA);
%[x0 x Lambda.lower Lambda.upper LB UB]

[results,table] = get_opf_results(om,x,Lambda,optns);
[results.et, results.f, results.success] = deal(et,f,success);

if optns.verify == 1
    success = verify_solutions(table,om,optns);
    fprintf(['\nSolutions verified: %i'],success);
end
%fd = fopen('output.txt','at');
%printpf(results,optns.fileID);'

% vv = get_idx(om);
% tab.Qg
% 180/pi*x(vv.i1.Va1:vv.iN.Va1)
table.Vm
table.Pg
table.Beta
table.Curtail
table.ExpCurtail
table.Qg
table.S
table.Slam
table.lamInfo

[h,g] = g_fcn(x);
g_dev = sum(abs(g))
h_dev = min(abs(h))

plot_results(table,optns);