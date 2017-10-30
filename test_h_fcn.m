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

%% generator options
optns.gen.fixPg = [1 3]; % generators which is fixed for stressed cases
optns.gen.fixQg = [];

optns.gen.maxPg = [1 3]; % generators for which to max production (must be fixed)
optns.gen.maxPgLim = [3000];

optns.branch.limit = 0; % turn on/off branch limits 
optns.branch.rateA = []; % branch limits

%optns.load.loadIncreaseArea = 1; % areas where to increase load
optns.bus.loadIncrease = [5]; % buses with load increase for contingencies
%% matpower options
optns.mpopt = mpoption();

optns.mpopt.pf.enforce_q_lims = 1;
optns.mpopt.opf.flow_lim = 'P';

%% solver options
foptions = optimoptions('fmincon','Algorithm','interior-point','GradObj','on','GradConstr','on');

foptions.Display = 'off'; % off, testing, iter
foptions.TolCon = 1e-10; % high value may give non-zero lagrange multipliers also for inactive constraints
foptions.TolFun = 1e2;
foptions.TolX = 1e-10;
foptions.MaxIter = 100;
%foptions.Algorithm = 'interior-point'; 

%optns.mpopt.opf.violation = 1e-3;
%% CHECK OPTIONS
optns = check_opf_options(optns);

%% SETUP MATPOWER CASE
mpc = setup_mpc(optns);

%mpc.gen(:,PG) = [0;0;0];

%% CONTINGENCIES
mpc = setup_contingencies(mpc,optns); 

%% convert to internal indexing
mpc = ext2int(mpc);


%% SETUPT OPTIMIZATION OBJECT
om = setup_opf(mpc,optns);

% build cost parameters, for Matpower hessian
om = build_cost_params(om);

% objective function
f_fcn = @(x)vscopf_f_maxPg(x, om);
% constraint function
g_fcn = @(x)vscopf_g(x, om, optns.mpopt);


[x0,LB,UB] = getv(om);

% solve using fmincon, without hessian
t0 = clock();
[x, f, success, Output, Lambda, grad, hessian_fmincon] = ...
  fmincon(f_fcn, x0, [], [], [], [], LB, UB, g_fcn, foptions);
et = etime(clock,t0);

% evaluate hessian function
hessian_fcn = vscopf_h(x,Lambda,om,optns.mpopt);
%function Lxx = opf_hessfcn(x, lambda, cost_mult, om, Ybus, Yf, Yt, mpopt, il)

max_diff = max(abs(hessian_fcn - hessian_fmincon)./ hessian_fmincon);

%Lxx_orig = opf_hessfcn(x,Lambda,0,om,Ybus,Yf,Yf,optns.mpopt,1:size(mpc.bus,1)); OK!