close all;
clear;
define_constants;
vscopf_define_constants;

%% OPTIONS
optns = struct();


optns.case = 'case9';
optns.contingencies = 'case9_contingencies';
optns.stabilityMargin = 0;

%% generator options
optns.gen.fixPg = [1]; % generators which is fixed for stressed cases
optns.gen.maxPg = [1]; % generators for which to max production (must be fixed)
optns.gen.maxPgLim = [3000];

optns.branch.limit = 0; % turn on/off branch limits 
optns.branch.rateA = []; % branch limits

optns.load.loadIncreaseArea = 1; % areas where to increase load
%% matpower options
optns.mpopt = mpoption();

optns.mpopt.pf.enforce_q_lims = 1;
optns.mpopt.opf.flow_lim = 'P';

%% solver options
foptions = optimoptions('fmincon');

foptions.Display = 'iter'; % off, testing, iter
foptions.GradObj = 'on';
foptions.GradConstr = 'on';
foptions.TolCon = 1e-10; % high value may give non-zero lagrange multipliers also for inactive constraints
foptions.TolFun = 1e2;
foptions.TolX = 1e-15;
foptions.MaxIter = 100;
%foptions.Algorithm = 'interior-point'; 

%optns.mpopt.opf.violation = 1e-3;
%% CHECK OPTIONS
check_opf_options(optns);

%% SETUP MATPOWER CASE
mpc = setup_mpc(optns);

% convert to internal indexing
mpc = ext2int(mpc);

%% CONTINGENCIES
mpc.contingencies = setup_contingencies(optns);

%% SETUPT OPTIMIZATION OBJECT
om = setup_opf(mpc,optns);

% objective function
f_fcn = @(x)vscopf_f_maxPg(x, om);
% constraint function
[Ybus,Yf,Yt] = makeYbus(mpc);
g_fcn = @(x)vscopf_g(x, om, Ybus, Yf, Yt, optns.mpopt);

[x0,LB,UB] = getv(om);

mpc = runpf(mpc,optns.mpopt);
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

%[x0 x Lambda.lower Lambda.upper LB UB]

[results,tab] = get_opf_results(om,x,Lambda);
[results.et, results.f, results.success] = deal(et,f,success);

%fd = fopen('output.txt','at');
printpf(results);