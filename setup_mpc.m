function mpc = setup_mpc(optns)

define_constants;

mpc = eval([optns.case]);



%% put extra generator information into mpc struct



% gen2
mpc.gen2 = [];
for i=1:size(mpc.gen,1)
    mpc.gen2(i,[GEN_BUS PMAXIMIZE PFIX]) = [mpc.gen(i,GEN_BUS) ismember(i,optns.gen.maxPg) ismember(i,optns.gen.fixPg)];
end

% change limits of generators which are being maximized
mpc.gen(mpc.gen2(:,PMAXIMIZE)==1,PMAX) = optns.gen.maxPgLim;


% bus2
mpc.bus2 = [];
for i=1:size(mpc.bus,1)
    if sum(mpc.bus(i,BUS_AREA) == optns.load.loadIncreaseArea) > 0
        mpc.bus2(i,LOAD_INCREASE_AREA) = 1;
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