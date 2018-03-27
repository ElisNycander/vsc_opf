function [mpc, rescase, restab, Output, om] = run_vscopf(mpc,optns,foptions,mpc0)
% Solve vscopf.
% mpc - actual case, used to calculate constraint functions (e.g. Pg base)
% mpc0 - initial values, used to setup optimization model object 

%% for script
% clear;
% load('randomize_initial_data.mat');
% optns.timer = 1;
% 
% foptions = optimoptions('fmincon','Algorithm','interior-point','GradObj','on','GradConstr','on');
% foptions.Display = 'off'; % off, testing, iter
% foptions.TolCon = 1e-10; % high value may give non-zero lagrange multipliers also for inactive constraints
% foptions.TolFun = 1e0;
% foptions.TolX = 1e-10;
% foptions.MaxIter = 2e3;
% foptions.MaxFunctionEvaluations = 5e3;

%% function
define_constants;
vscopf_define_constants;

%% RANDOMIZE INITIAL POINT

%% RUN INITIAL POWER FLOW FOR BASE CASE

%% matpower options
%optns.mpopt = mpoption();
%optns.mpopt.opf.flow_lim = 'S';

% do pfs quietly
%optns.mpopt.out.all = 0;
%optns.mpopt.verbose = 0;
%optns.mpopt.pf.enforce_q_lims = 0;
mpci = runpf(mpc,optns.mpopt);

%mpc = mpci; % may change bus types, including slack
%mpci = runpf(mpc); % don't enforce q limits
if optns.useInitialPF
    assert(mpci.success == 1,'Initial power flow was unsuccessful');
    % copy base case values, for initial guess
    mpc.gen(:,[QG PG]) = mpci.gen(:,[QG PG]);
    %mpc.gen(:,PG) = mpci.gen(:,PG);
    mpc.bus(:,[VM VA]) = mpci.bus(:,[VM VA]);
    assert(min(mpc.gen(:,PG) >= 0),'At least one generator with negative Pg for base case')
    %mpc.bus(:,BUS_TYPE) = mpci.bus(:,BUS_TYPE);
end
%save('mpc_basecase.mat','mpci');

%% convert to internal indexing
mpc = ext2int(mpc);
if nargin > 3
mpc0 = ext2int(mpc0);
end
%% CONTINGENCIES
mpc = setup_contingencies(mpc,optns); 

%mpc0.contingencies = mpc.contingencies;
%% SETUP OPTIMIZATION OBJECT
if nargin > 3
    om = setup_opf_x0(mpc,optns,mpc0);
else
    om = setup_opf(mpc,optns);
end
    
    
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

t0 = clock();
[x, f, success, Output, Lambda] = ...
    fmincon(f_fcn, x0, A, U, [], [], LB, UB, g_fcn, foptions);
et = etime(clock,t0);
if optns.timer == 1
    fprintf(['Time: %0.3f s'],et)
    fprintf(['\nObjective function: %0.1f\n'],f*mpc.baseMVA);
end
%[x0 x Lambda.lower Lambda.upper LB UB]

[h,g] = g_fcn(x);
g_dev = sum(abs(g));
h_dev = max(h);

% collect results in tables
[rescase,restab] = get_opf_results(om,x,Lambda,optns);

% verify solutions
restab.solutions_verified = verify_solutions_sbus(restab,om,optns);
restab.Output = Output;
[restab.et, restab.f, restab.success, restab.g_dev, restab.h_dev] =  ...
    deal(et,f,success, g_dev, h_dev);



