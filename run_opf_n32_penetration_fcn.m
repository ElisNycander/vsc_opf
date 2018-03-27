function exitflag = run_opf_n32_penetration_fcn(caseName,windScenario,usePQConstraints,enableQWind,penetration)

%close all;
%clear;
define_constants;
vscopf_define_constants;

%penetration = 0.5:0.5:0.9

%% OPTIONS
optns = struct();

[optns.caseName, optns.windScenario, optns.gen.usePQConstraints, optns.QWind, optns.penetrationLevel] = ...
   deal(caseName,windScenario, usePQConstraints, enableQWind, penetration);

% optns.caseName = 'case_default';
% optns.windScenario = 'D1'; % 
% optns.gen.usePQConstraints = 0;
% optns.QWind = 0;
% optns.penetrationLevel = 0.3;

%%

optns.outputFile = '';

optns.caseFile = 'Nordic32';
%eval(['mpc=' optns.caseFile ';']); % get case
mpc = Nordic32;

optns.contingencyFile = 'n32_contingencies';
%optns.stabilityMargin = 0.1;

optns.hessian = 0;
optns.verify = 1; % verify_solutions
optns.verifyPQconversion = 0; % convert PV buses with generators at Q-limits to PQ in power flow

optns.lambdaTolerance = 1e-2; % round smaller lambda to 0 for tables

optns.useInitialPF = 1; % initial power flow must be solvable

optns.useOlaussonScenarios = 1; % Scenario from Olausson (2015)
optns.setPenetration = 1; % Manually set penetration ratio of system for base case to this value

optns.replaceGeneration = 1; % replace wind with synchronous generation in base case
optns.generationReplacementTol = 5; % in MW

optns.lowLoadScenario = 0; % reduce load for low-load scenario

highLoad = 22e3; % high load (in MW)
lowLoad = 10e3; % low load (in MW)

slackFactor = 1; % scale down generation, to get positive production at slack bus in base case

optns.gen.optimizeBaseP = 0; 
optns.gen.fixBaseWind = 1; % fix curtailable P for base case (only when optimizBbaseP) - NOT IMPLEMENTED


optns.saveFigures = 1;
optns.saveData = 1;

% activate/deactivate given wind power scenarios, if active the scenarios
% will be constructed as [wind scenarios x contingencies]
optns.gen.useWindScenarios = 1;

optns.useWindVariance = 1;
optns.gen.windVariance = 0.0083;
nr_std = 1.96;

optns.transferCorridors = {
    [4011 4071;4012 4071], [4031 4041;4032 4044;4032 4042;4021 4042]
}; % transfer over these lines will be summed (active power flow P)   
optns.externalBuses = [4071 4072]; % used when computing wind pentration in plot_results
%%
mpc = n32_define_areas(mpc);

% run load flow for original system
op = mpoption();
op.pf.enforce_q_lims = 0;
op.verbose = 0;
op.out.all = 0;
mpco = runpf(mpc,op);


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
%     1011   0 0   50      -50       1   100     1       1e3     0
%     2031   0 0   50      -50       1   100     1       1e3     0
%     41     0 0   50      -50       1   100     1       1e3     0
%     1041   0 0   50      -50       1   100     1       1e3     0  
  
% choose some nodes
% North
    4011   0 0   0      -0       1   100     1       1e3     0
    1011   0 0   0      -0       1   100     1       1e3     0
    1014   0 0   0      -0       1   100     1       1e3     0
    1012   0 0   0      -0       1   100     1       1e3     0
    4022   0 0   0      -0       1   100     1       1e3     0
    4021   0 0   0      -0       1   100     1       1e3     0
    2031   0 0   0      -0       1   100     1       1e3     0
% Central and Southwest
    1041   0 0   0      -0       1   100     1       1e3     0  
    1045   0 0   0      -0       1   100     1       1e3     0  
    4051  0 0   0      -0       1   100     1       1e3     0  
    4046  0 0   0      -0       1   100     1       1e3     0  
    4061  0 0   0      -0       1   100     1       1e3     0  
    1044  0 0   0      -0       1   100     1       1e3     0  
    4062  0 0   0      -0       1   100     1       1e3     0  
    
    % choose some nodes, ordered by increasing bus idx
% North
%     1011   0 0   0      -0       1   100     1       1e3     0
%     1012   0 0   0      -0       1   100     1       1e3     0
%     1014   0 0   0      -0       1   100     1       1e3     0
%     1041   0 0   0      -0       1   100     1       1e3     0
%     1044   0 0   0      -0       1   100     1       1e3     0
%     1045   0 0   0      -0       1   100     1       1e3     0
%     2031   0 0   0      -0       1   100     1       1e3     0
% % Central and Southwest
%     4011   0 0   0      -0       1   100     1       1e3     0  
%     4021   0 0   0      -0       1   100     1       1e3     0  
%     4022  0 0   0      -0       1   100     1       1e3     0  
%     4046  0 0   0      -0       1   100     1       1e3     0  
%     4051  0 0   0      -0       1   100     1       1e3     0  
%     4061  0 0   0      -0       1   100     1       1e3     0  
%     4062  0 0   0      -0       1   100     1       1e3     0  
];

if optns.lowLoadScenario && ~ optns.setPenetration
    optns.gen.extra(:,QMIN) = -100;
    optns.gen.extra(:,QMAX) = 100;
end

% Sweden peak load ~ 25000 MW, min load ~ 9000 MW
% subtract 2300 MW from external buses
optns.gen.scaleFactor = ( sum(mpc.bus(:,PD)) ...
        - sum( mpc.bus(mpc.bus(:,BUS_AREA) == 4,PD)) ) / highLoad; 


% for PQ capability
%optns.gen.pqFactor = ones(size(mpc.gen,1)+size(optns.gen.extra,1),1); % constraint omitted for 0 value
%optns.gen.pqFactor = zeros(size(mpc.gen,1)+size(optns.gen.extra,1),1);
% optns.gen.pqFactor = [
% %     ones(3,1)
% %     zeros(10,1)
% %     1
% %     1
% %     0.01
% %     2
% %     zeros(10,1)
% zeros(13,1)
% 100
% zeros(13,1)
%    % zeros(24,1)
%    % 2*ones(23,1)
%    % zeros(4,1)
% ];
if ~optns.QWind
    optns.gen.pqFactor = [
        1*ones(12,1)
        0 % note gen 13 is synchronous condenser
        1*ones(10,1)
        0*ones(14,1)
        ];
else
    optns.gen.pqFactor = [
        1*ones(12,1)
        0 % note gen 13 is synchronous condenser
        1*ones(10,1)
        0*ones(14,1)
        %sqrt(1-0.9^2)/0.9*ones(14,1) % PF 0.9
        ];
end
%mpc.gen(11,PG) = mpc.gen(11,PG)-50;
%mpc.gen(10,PG) = mpc.gen(10,PG)-50;
%optns.gen.compBuses = [4021 4012]; % generators at these buses will have their PG decreased to compensate for extra generation in base case
%%
gen_bus_2_digits = rem(mpc.gen(:,GEN_BUS),100);
fixedP_boolean = and(gen_bus_2_digits > 40, gen_bus_2_digits < 70);

%% the active power of generators can either be:
% 1 Variable - can vary freely for all contingencies
% 2 Fixed - set by base case power flow
% 3 Curtailable - can be curtailed relative to base case in contingencies
% Note: Variable is default
%optns.gen.fixedP = [15 16 17 18];
%optns.gen.fixedP = [15:18 20 21]; % all nuclear reactors

% only vary generators in North
optns.gen.fixedP = find(gen_bus_2_digits > 40);

%optns.gen.fixedP = find(fixedP_boolean);
optns.gen.curtailableP = [24:24+size(optns.gen.extra,1)-1];
optns.gen.variableP = [];
% include non-variable p in optimization or not, may take given values as "market outcome"


% wind power scenarios - stored in columns, with as many rows as there are
% curtailable generators, or just one row, in which case the same scenario
% is applied to all curtailable generators
% optns.gen.windScenarios = [
%     2000
% ];

if optns.useOlaussonScenarios
    % use Olausson scenarios
    optns.gen.capacityTable = olausson_scenarios(optns.windScenario);
    
    % find area of new generators
    optns.gen.genAreaExtra = zeros(size(optns.gen.extra,1),1);
    for i=1:length(optns.gen.genAreaExtra)
        bus = optns.gen.extra(i,GEN_BUS);
        optns.gen.genAreaExtra(i) = mpc.bus(mpc.bus(:,BUS_I)==bus,BUS_AREA);
    end
    % find area of old generators
    optns.gen.genArea = zeros(size(mpc.gen,1),1);
    for i=1:length(optns.gen.genArea)
        bus = mpc.gen(i,GEN_BUS);
        optns.gen.genArea(i) = mpc.bus(mpc.bus(:,BUS_I)==bus,BUS_AREA);
    end
    
    % preallocate
    optns.gen.windScenarios = zeros(size(optns.gen.extra,1),1);
    
    % divide generators into groups
    northIdxExtra = optns.gen.genAreaExtra == 1;
    southIdxExtra = or(optns.gen.genAreaExtra == 2,optns.gen.genAreaExtra == 3);
    northIdxSynch = optns.gen.genArea == 1;
    southIdxSynch = or(optns.gen.genArea == 2, optns.gen.genArea == 3);
    
    if ~ optns.setPenetration % penetration level from case is taken
    % SE1, SE2 -> North
    optns.gen.windScenarios(northIdxExtra) = optns.gen.scaleFactor * ...
        sum(optns.gen.capacityTable(1:2))/sum(northIdxExtra);
    % SE3, SE4 -> Central and South
    optns.gen.windScenarios(southIdxExtra) = optns.gen.scaleFactor * ...
        sum(optns.gen.capacityTable(3:4))/sum(southIdxExtra);
    
    if optns.lowLoadScenario % create low-load scenario
        
        % scale down load (not external)
        mpc.bus(mpc.bus(:,BUS_AREA) ~= 4,[PD QD]) = mpc.bus(mpc.bus(:,BUS_AREA)~= 4,[PD QD])*lowLoad/highLoad;
        
        % deactivate some nuclear
        deactivatedGens = [18 21 16];
        mpc.gen(deactivatedGens,[PG QG PMIN PMAX QMIN QMAX]) = 0;
        
        
        % put base scenario into original PF
        optns.gen.extra(:,PG) = optns.gen.windScenarios(:,1);
        
        % find conventional generators (internal)
        idxConvIntGen = gen_idx( mpc.bus(mpc.bus(:,BUS_AREA)~=4,BUS_I),mpc);
        
        % surplus = wind + conventional (internal) - load (internal)
        windTot = sum(optns.gen.extra(:,PG));
        loadTot = sum(mpc.bus(mpc.bus(:,BUS_AREA)~=4,PD));
        %surplus = sum(optns.gen.extra(:,PG)) + sum(mpc.gen(idxConvIntGen,PG)) - ...
        %    sum(mpc.bus(mpc.bus(:,BUS_AREA) ~= 4,PD));
        
        % total wind must be less than interal load
        assert(windTot < loadTot,...
            'Total wind must be less than internal load');
        
        % scale down other generation
        mpc.gen(idxConvIntGen,[PG QG]) = mpc.gen(idxConvIntGen,[PG QG]) * ...
                                          (loadTot - windTot) / sum(mpc.gen(idxConvIntGen,PG));
        
        % how much generation is needed to match load
    end
    
    
    else % set penetration level
        
        idxSweGen = gen_bus_2_digits < 70;
        totGen = sum(mpc.gen(idxSweGen,PG));
        totWind = optns.penetrationLevel * totGen;
        
        nshare = sum(optns.gen.capacityTable(1:2))/optns.gen.capacityTable(5);
        sshare = sum(optns.gen.capacityTable(3:4))/optns.gen.capacityTable(5);
        
        optns.gen.windScenarios(northIdxExtra) = totWind * nshare / sum(northIdxExtra);
        optns.gen.windScenarios(southIdxExtra) = totWind * sshare / sum(southIdxExtra);
        
        
    end
    
    
    if optns.replaceGeneration
        
        if optns.useWindVariance
            % use 95% confidence interval for wind production
            optns.gen.extra(:,PG) = optns.gen.windScenarios(:,1) * ...
                                (1-sqrt(optns.gen.windVariance)*nr_std);
        else
            % put base scenario into original PF
            optns.gen.extra(:,PG) = optns.gen.windScenarios(:,1);
        end
        
        % reduce other generation
        northWindGen = sum(northIdxExtra .* optns.gen.extra(:,PG));
        southWindGen = sum(southIdxExtra .* optns.gen.extra(:,PG));
        
        northGen = sum(northIdxSynch .* mpco.gen(:,PG));
        southGen = sum(southIdxSynch .* mpco.gen(:,PG));
        
        
        
        
        % check if wind power in south replaces all conventional generation
        if southWindGen > southGen - optns.generationReplacementTol
            northWindGen = northWindGen + southWindGen - southGen;
            % set all generation in south to 0
            mpc.gen(southIdxSynch,PG) = 0;
        else
            mpc.gen(southIdxSynch,PG) = mpc.gen(southIdxSynch,PG).* (1 - southWindGen/southGen);
        end
        mpc.gen(northIdxSynch,PG) = mpc.gen(northIdxSynch,PG).* (1 - northWindGen/northGen);
        
        
        mpc.gen(:,PG) = mpc.gen(:,PG)*slackFactor;
    end
else
    optns.gen.windScenarios = [500];
end

% increase and decrease in wind production
if optns.useWindVariance
    optns.gen.windScenarios = [
        optns.gen.windScenarios * [1 1-2*nr_std*sqrt(optns.gen.windVariance)]
        ];
end

newScenario = zeros(size(optns.gen.windScenarios(:,1)));
newScenario(northIdxExtra) = optns.gen.windScenarios(northIdxExtra,1);
newScenario(southIdxExtra) = (1-2*nr_std*sqrt(optns.gen.windVariance))*...
                            optns.gen.windScenarios(southIdxExtra,1);
% make new scenario: increase in north, decrease in south
optns.gen.windScenarios = [ optns.gen.windScenarios newScenario ];

%northSynchGen = sum(mpc.gen

% probabilities for wind scenarios, empty means all scenarios equally
% likely
optns.gen.windProbabilities = [
];



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
rateA(rateA == 4000) = 1000;
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
optns.branch.rateA = rateA;


optns.branch.duplicate = []; % duplicate these branches

%% load increase
%optns.load.loadIncreaseArea = 1; % areas where to increase load
optns.bus.loadIncrease = []; % buses with load increase for contingencies



%% solver options
foptions = optimoptions('fmincon','Algorithm','interior-point','GradObj','on','GradConstr','on');

foptions.Display = 'off'; % off, testing, iter
foptions.TolCon = 1e-10; % high value may give non-zero lagrange multipliers also for inactive constraints
foptions.TolFun = 1e0;
foptions.TolX = 1e-10;
foptions.MaxIter = 2e3;
foptions.MaxFunctionEvaluations = 10e3;

%% CHECK OPTIONS
optns = check_opf_options(mpc,optns);

%% SETUP MATPOWER CASE
mpc = setup_mpc(mpc,optns);

% set min and max values for PG (not present in N32)
mpc.gen(:,[PMAX PMIN]) = [mpc.gen(:,MBASE) zeros(size(mpc.gen,1),1)];


if optns.lowLoadScenario && ~optns.setPenetration   
    % change Q-limits
    mpc.gen(:,QMIN) = -1e4;
    % set limits for inactive generators
    mpc.gen(deactivatedGens,[PG PMIN PMAX QMIN QMAX]) = 0;
end

if optns.QWind
    %    mpc.gen(optns.gen.curtailableP,QMIN) = - mpc.gen(optns.gen.curtailableP,PG)* sqrt(1-0.81)/0.9;
    %    mpc.gen(optns.gen.curtailableP,QMAX) = + mpc.gen(optns.gen.curtailableP,PG)* sqrt(1-0.81)/0.9;
    
    mpc.gen(optns.gen.curtailableP,QMIN) = - optns.gen.windScenarios(:,1) * sqrt(1-0.81)/0.9;
    mpc.gen(optns.gen.curtailableP,QMAX) =   optns.gen.windScenarios(:,1) * sqrt(1-0.81)/0.9;
end

% fix P limits for synchronous condenser
mpc.gen(find(mpc.gen(:,GEN_BUS)==4041),PMAX) = 0;



%% RUN INITIAL POWER FLOW FOR BASE CASE

%% matpower options
optns.mpopt = mpoption();

optns.mpopt.opf.flow_lim = 'S';

% do pfs quietly
optns.mpopt.out.all = 0;
optns.mpopt.verbose = 4;
optns.mpopt.pf.enforce_q_lims = 0;
mpci = runpf(mpc,optns.mpopt);

%mpc = mpci; % may change bus types, including slack
%mpci = runpf(mpc); % don't enforce q limits
if optns.useInitialPF
    assert(mpci.success == 1,'Initial power flow was unsuccessful');
    % copy base case values, for initial guess
    mpc.gen(:,[QG PG]) = mpci.gen(:,[QG PG]);
    mpc.bus(:,[VM VA]) = mpci.bus(:,[VM VA]);
    assert(min(mpc.gen(:,PG) >= 0),'At least one generator with negative Pg for base case')
    %mpc.bus(:,BUS_TYPE) = mpci.bus(:,BUS_TYPE);
end
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

[rescase,restab] = get_opf_results(om,x,Lambda,optns);
[rescase.et, rescase.f, rescase.success] = deal(et,f,success);

if optns.verify == 1
    check_solutions = verify_solutions(restab,om,optns);
    fprintf(['\nSolutions verified: %i'],sum(check_solutions)==mpc.contingencies.N);
end


%fd = fopen('output.txt','at');
%printpf(rescase,optns.fileID);'

% vv = get_idx(om);
% tab.Qg
% 180/pi*x(vv.i1.Va1:vv.iN.Va1)
%restab.Vm
restab.Pg
%restab.Beta
%restab.Curtail
%restab.ExpCurtail
restab.Qg
restab.Flow
%restab.Flowlam
%restab.lamInfo
restab.PQ
restab.Pflow
restab.Qflow
restab.transferFrom
restab.transferTo

[h,g] = g_fcn(x);
g_dev = sum(abs(g))
h_dev = max(h)

%plot_results(restab,optns);
%close all;
if optns.saveData == 1
    save(['data/' optns.caseName '.mat'],'rescase','restab','om','x','Lambda','optns','Output','f');
end


exitflag = success;
