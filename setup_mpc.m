function mpc = setup_mpc(optns)

define_constants;

mpc = eval([optns.caseFile]);

%% extra generators
ng = size(mpc.gen,1);
ng_extra = size(optns.gen.extra,1);
nrow = size(optns.gen.extra,2);
% add extra generators
mpc.gen(ng+1:ng+ng_extra,1:nrow) = optns.gen.extra;
if ~isfield(optns.gen,'extra_cost') 
   mpc.gencost(ng+1:ng+ng_extra,1:size(mpc.gencost,2)) = 0; 
else
   mpc.gencost(ng+1:ng+ng_extra,1:size(optns.gen.extra_cost,2)) = optns.gen.extra_cost;
end

%% extra lines


%% put extra generator information into mpc struct



% gen2
mpc.gen2 = [];
for i=1:size(mpc.gen,1)
    mpc.gen2(i,[GEN_BUS PMAXIMIZE]) = ...
        [mpc.gen(i,GEN_BUS) ...
         ismember(i,optns.gen.maxPg) ...
      ];
    if ismember(i,optns.gen.fixedP)
        mpc.gen2(i,PTYPE) = PFIX;
    elseif ismember(i,optns.gen.curtailableP)
        mpc.gen2(i,PTYPE) = PCUR;
    else
        mpc.gen2(i,PTYPE) = PVAR;
    end
end

% change limits of generators which are being maximized
mpc.gen(mpc.gen2(:,PMAXIMIZE)==1,PMAX) = optns.gen.maxPgLim;


% bus2
mpc.bus2 = [];
for i=1:size(mpc.bus,1)
    if ismember(i,optns.bus.loadIncrease)
%    if sum(mpc.bus(i,BUS_AREA) == optns.load.loadIncreaseArea) > 0
        mpc.bus2(i,LOAD_INCREASE_AREA) = 1;
    else
        mpc.bus2(i,LOAD_INCREASE_AREA) = 0;
    end
end


mpc.stabilityMargin = optns.stabilityMargin;


% branch flow limits
if ~optns.branch.limit
    mpc.branch(:,RATE_A) = 0;
else
    if ~isempty(optns.branch.rateA)
        mpc.branch(:,RATE_A) = optns.branch.rateA;
    end
end

%% duplicate branches

for i=1:size(optns.branch.duplicate)
    mpc.branch = [mpc.branch; mpc.branch(optns.branch.duplicate(i),:)];
end
