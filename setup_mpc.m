function mpc = setup_mpc(mpc,optns)
% Adds the following matrices to the mpc struct:
% gen2
% bus2
%
% Also modifies branch matrix
% Note: Uses external indexing

define_constants;

%mpc = eval([optns.caseFile]);

%% extra generators
ng = size(mpc.gen,1);
ng_extra = size(optns.gen.extra,1);
nrow = size(optns.gen.extra,2);
% add extra generators
mpc.gen(ng+1:ng+ng_extra,1:nrow) = optns.gen.extra;

% only do this if cost matrix exists
if isfield(mpc,'gencost')
    if ~isfield(optns.gen,'extra_cost')
        % add zero matrix
        mpc.gencost(ng+1:ng+ng_extra,1:size(mpc.gencost,2)) = 0;
    else
        % add new cost information
        mpc.gencost(ng+1:ng+ng_extra,1:size(optns.gen.extra_cost,2)) = optns.gen.extra_cost;
    end
end

% define all generators as active 
inactiveGen = find(mpc.gen(:,GEN_STATUS) == 0);
if ~isempty(inactiveGen)
    warning('Inactive generators were set to active in base scenario')
    mpc.gen(inactiveGen,GEN_STATUS) = 1;
end

% set VG for new generators to VG for old generators at that bus
for i=1:ng_extra
   thisGenIdx = ng+i;
   thisBus = mpc.gen(thisGenIdx,GEN_BUS);
   sameBusGenIdx = find(mpc.gen(1:ng,GEN_BUS)==thisBus);
   if ~isempty(sameBusGenIdx)
        mpc.gen(thisGenIdx,VG) = mpc.gen(sameBusGenIdx(1),VG); % choose first generator at this bus
        warning('VG set point at new generators changed to agree with old generators')
   end
end

% set all buses with generators to PV buses
for i=1:size(mpc.bus,1)
    if ~ isempty( find( mpc.gen(:,GEN_BUS) == mpc.bus(i,BUS_I),1 ) )
        if mpc.bus(i,BUS_TYPE) == PQ % Change PQ buse to PV
            mpc.bus( i,BUS_TYPE ) = PV;
        end
    end
end
    
% % compensate PG to maintain active power balance
% % find all generators at buses where to perform compensation
% compGen = [];
% for i=1:length(optns.compBuses)
%     compGen = [compGen find(mpc.gen(:,GEN_BUS)==optns.compBuses(i))];
% end
% if isempty(compGen)
%    % compensate att  
%     
% end


%% extra lines


%% put extra generator information into mpc struct



% gen2
mpc.gen2 = [];
for i=1:size(mpc.gen,1)
    mpc.gen2(i,[GEN_BUS]) = ...
        [mpc.gen(i,GEN_BUS) ...
      ];
    if ismember(i,optns.gen.fixedP)
        mpc.gen2(i,PTYPE) = PFIX;
    elseif ismember(i,optns.gen.curtailableP)
        mpc.gen2(i,PTYPE) = PCUR;
    else
        mpc.gen2(i,PTYPE) = PVAR;
    end
end

if ~isempty(optns.gen.pqFactor)
    mpc.gen2(:,PQ_FACTOR) = optns.gen.pqFactor;
else
    mpc.gen2(:,PQ_FACTOR) = 0;
end

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


%mpc.stabilityMargin = optns.stabilityMargin;


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
