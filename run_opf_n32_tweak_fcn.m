function exitflag = run_opf_n32_tweak_fcn(caseName,windScenario,usePQConstraints,enableQWind,penetration,plotResults)
optns = struct();

[optns.caseName,optns.windScenario, optns.gen.usePQConstraints, optns.QWind, optns.penetrationLevel] = ...
   deal(caseName,windScenario, usePQConstraints, enableQWind, penetration);
%% for script
% clear;
% close all;
% 
% optns.caseName = 'E1_default';
% optns.windScenario = 'E1';
% 
% 
% optns.gen.usePQConstraints = 1;
% optns.QWind = 0;
% 
% plotResults = 1;
% 
% optns.penetrationLevel = 0.5;


%% function
define_constants;
vscopf_define_constants;
 optns.verifyTol = 1; 
optns.timer = 1;
optns.baseCaseTransferIncrease = 400;
optns.reduceAllGeneration = 0; % reduce all generation to compensate for wind power, i.e. i
                                % increased transfers
%% OPTIONS
optns.mpopt = mpoption();
optns.mpopt.opf.flow_lim = 'S';

% do pfs quietly
optns.mpopt.out.all = 1;
optns.mpopt.verbose = 2;
optns.mpopt.pf.enforce_q_lims = 0;


optns.outputFile = '';

optns.caseFile = 'Nordic32';
mpc = Nordic32;

optns.contingencyFile = 'n32_contingencies';
%optns.stabilityMargin = 0.1;

optns.hessian = 0;
optns.verify = 1; % verify_solutions
optns.verifyPQconversion = 0; % convert PV buses with generators at Q-limits to PQ in power flow

optns.lambdaTolerance = 1e-3; % round smaller lambda to 0 for tables

optns.useInitialPF = 1; % initial power flow must be solvable


%optns.windArea = 2; 
% 1 - North
% 2 - North + Central

%optns.penetrationLevel = 0.7;
%optns.replaceGeneration = 1; % replace wind with synchronous generation in base case
%optns.generationReplacementTol = 5; % in MW


highLoad = 22e3; % high load (in MW)
lowLoad = 10e3; % low load (in MW)

slackFactor = 1; % scale down generation, to get positive production at slack bus in base case

optns.gen.optimizeBaseP = 0; 
optns.gen.fixBaseWind = 1; % fix curtailable P for base case (only when optimizBbaseP) - NOT IMPLEMENTED

optns.saveFigures = 0;
optns.saveData = 1;

% activate/deactivate given wind power scenarios, if active the scenarios
% will be constructed as [wind scenarios x contingencies]
optns.gen.useWindScenarios = 1;

optns.useWindVariance = 1;
optns.gen.windVariance = 0.0083;
nr_std = 1.96; % for 95% confidence interval

optns.transferCorridors = {
    [4031 4041;4032 4044;4032 4042;4021 4042], [4011 4071;4012 4071]
}; % transfer over these lines will be summed (active power flow P)   
optns.externalBuses = [4071 4072]; % used when computing wind pentration in plot_results

optns.branch.limit = 0; % turn on/off branch limits % Note: if optimizeBaseP = 0 limits are only added for contingencies
optns.branch.limit_base_case = 0;

%% solver options
foptions = optimoptions('fmincon','Algorithm','interior-point','GradObj','on','GradConstr','on');

foptions.Display = 'off'; % off, testing, iter
foptions.TolCon = 1e-10; % high value may give non-zero lagrange multipliers also for inactive constraints
foptions.TolFun = 1e0;
foptions.TolX = 1e-10;
foptions.MaxIter = 2e3;
foptions.MaxFunctionEvaluations = 5e3;

%%

mpc = n32_define_areas(mpc);

% run load flow for original system
op = mpoption();
op.pf.enforce_q_lims = 0;
op.verbose = 0;
op.out.all = 0;
mpco = runpf(mpc,op);

% set slack PG = 0
mpc.gen( mpc.gen(:,GEN_BUS) == mpc.bus( mpc.bus(:,BUS_TYPE) == 3,BUS_I ) ,PG) = 0;

totalGenOriginal = sum(mpco.gen(:,PG));

%% generator options
% extra generators           
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	 

optns.gen.extra = wind_scenario(optns);
nWindFarms = size(optns.gen.extra,1);

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

% divide generators into groups
northIdxExtra = optns.gen.genAreaExtra == 1;
southIdxExtra = or(optns.gen.genAreaExtra == 2,optns.gen.genAreaExtra == 3);
northIdxSynch = optns.gen.genArea == 1;
southIdxSynch = or(optns.gen.genArea == 2, optns.gen.genArea == 3);

ncIdxSynch = or(optns.gen.genArea == 1, optns.gen.genArea == 2);

%% increase transfers in base case

genP = optns.gen.genArea == 1;
genN = optns.gen.genArea == 2;

mpc.gen( genP,PG ) = mpc.gen( genP,PG ) * ( ...
            1 + optns.baseCaseTransferIncrease / sum( mpc.gen ( genP,PG ) ) );
mpc.gen( genN,PG ) = mpc.gen( genN,PG ) * ( ...
            1 - optns.baseCaseTransferIncrease / sum( mpc.gen ( genN,PG ) ) );


% Sweden peak load ~ 25000 MW, min load ~ 9000 MW
% subtract 2300 MW from external buses
% optns.gen.scaleFactor = ( sum(mpc.bus(:,PD)) ...
%         - sum( mpc.bus(mpc.bus(:,BUS_AREA) == 4,PD)) ) / highLoad;


% for PQ capability

if optns.QWind
    optns.gen.pqFactor = [
        1*ones(12,1)
        0 % note gen 13 is synchronous condenser
        1*ones(10,1)
        sqrt(1-0.9^2)/0.9*ones(size(optns.gen.extra,1),1)
        ];
else
    optns.gen.pqFactor = [
        1*ones(12,1)
        0 % note gen 13 is synchronous condenser
        1*ones(10,1)
        0*ones(size(optns.gen.extra,1),1)
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
optns.gen.fixedP = [15:18 20 21 22 23]; % all nuclear reactors, and 4072
%optns.gen.fixedP = [];
% only vary generators in North
%optns.gen.fixedP = find(gen_bus_2_digits > 40);
%optns.gen.fixedP = find(gen_bus_2_digits > 70);

%optns.gen.fixedP = find(fixedP_boolean);
optns.gen.curtailableP = [24:24+nWindFarms-1];
optns.gen.variableP = [];
% include non-variable p in optimization or not, may take given values as "market outcome"


% wind power scenarios - stored in columns, with as many rows as there are
% curtailable generators, or just one row, in which case the same scenario
% is applied to all curtailable generators

% preallocate
optns.gen.windScenarios = zeros(size(optns.gen.extra,1),1);



idxSweGen = gen_bus_2_digits < 70;
%totGen = sum(mpc.gen(idxSweGen,PG)); % total synchronous generation
switch optns.windScenario
    
    case 'E1' % Only wind in North
        
%         if optns.reduceAllGeneration
%             idxReduceGen = ncIdxSynch;
%         else
%             idxReduceGen = northIdxSynch;
%         end
        northGen = sum( mpc.gen( northIdxSynch,PG ) );
        totGen = sum(mpc.gen( ncIdxSynch ,PG));
        totWind = optns.penetrationLevel * northGen; % total wind
        
        % distribute generation equally over wind farms
        nshare = sum(northIdxExtra)/length(northIdxExtra);
        sshare = sum(southIdxExtra)/length(southIdxExtra);
        
        optns.gen.windScenarios(northIdxExtra) = totWind * nshare / sum(northIdxExtra);
        optns.gen.windScenarios(southIdxExtra) = totWind * sshare / sum(southIdxExtra);
        
        % put wind into gen matrix
        optns.gen.extra(:,PG) = optns.gen.windScenarios(:,1);
        
        % reduce other generation
        if optns.reduceAllGeneration
            mpc.gen( ncIdxSynch,PG ) = (1-totWind/totGen) * mpc.gen( ncIdxSynch,PG );
        else
            mpc.gen( northIdxSynch,PG ) = (1-totWind/northGen) * mpc.gen( northIdxSynch,PG );
        end
        %mpc.gen(idxReduceGen,PG) = mpc.gen(idxReduceGen,PG)*slackFactor;
        
        
    case 'E2'
        % penetration based central and north
        totGenCN = sum(mpc.gen( ncIdxSynch , PG ));
        totWind = optns.penetrationLevel * totGenCN;
        
        % distribute generation equally over wind farms
        optns.gen.windScenarios(:) = totWind / nWindFarms;
        
        % put wind into gen matrix
        optns.gen.extra(:,PG) = optns.gen.windScenarios(:,1);
        
        % reduce other generation
        mpc.gen(ncIdxSynch,PG) = (1-totWind/totGenCN) * mpc.gen(ncIdxSynch,PG);
        
end

% increase and decrease in wind production
if optns.useWindVariance
    optns.gen.windScenarios = [
        optns.gen.windScenarios * [1-nr_std*sqrt(optns.gen.windVariance) 1+nr_std*sqrt(optns.gen.windVariance)]
        ];
end


% probabilities for wind scenarios, empty means all scenarios equally
% likely
optns.gen.windProbabilities = [
    ];



%% branch options

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
lines = optns.transferCorridors{1};
idxCritLines = line_idx(lines(:,1),lines(:,2),mpco);

rateA(rateA == 1000) = 350;
rateA(rateA == 2000) = 500;
%rateA(rateA == 4000) = 1400/2;
rateA(rateA == 4000) = 1400;
%rateA(idxCritLines) = 850;
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



%% CHECK OPTIONS
optns = check_opf_options(mpc,optns);

%% SETUP MATPOWER CASE
mpc = setup_mpc(mpc,optns);

% set min and max values for PG (not present in N32)
mpc.gen(:,[PMAX PMIN]) = [mpc.gen(:,MBASE) zeros(size(mpc.gen,1),1)];

if optns.QWind
    mpc.gen(optns.gen.curtailableP,QMIN) = - mpc.gen(optns.gen.curtailableP,PG)* sqrt(1-0.81)/0.9;
    mpc.gen(optns.gen.curtailableP,QMAX) = + mpc.gen(optns.gen.curtailableP,PG)* sqrt(1-0.81)/0.9;
end

% fix P limits for synchronous condenser
mpc.gen(find(mpc.gen(:,GEN_BUS)==4041),PMAX) = 0;


%% RUN VSCOPF

[mpc, rescase, restab, Output, om] = run_vscopf(mpc,optns,foptions);

%% Verify solutions
if optns.verify == 1
    success = restab.solutions_verified;
    fprintf(['\nSolutions verified: %i\n'],sum(success)==mpc.contingencies.N);
end

%% Plots results
if plotResults
    plot_results(restab,optns);
end

if optns.saveData == 1
    save(['data/' optns.caseName '.mat'],'Output','rescase','restab','om','optns');
end
exitflag = restab.success;