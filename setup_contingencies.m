function mpc = setup_contingencies(mpc,optns)
% Construct lists of active generators, active lines, and stacked
% Y-matrices for each scenario, to be used when evaluating constraint and 
% objective functions.
% Note: mpc is converted to internal ordering, but contingencies and
% settings are given in external ordering. 

define_constants;

list = eval(optns.contingencyFile);

if isempty(list) % minimal list if there are no contigencies
    list = zeros(1,4);  
end
N = size(list,1);
nContingencies = N;
% add column with equal probabilities for all contingencies (except base
% case)
if size(list,2) < CONT_PROB
    list = [list 1/N*ones(N,1)];
end

% [contingencies x wind scenarios]
nScenarios = size(optns.gen.windScenarios,2);

%% create ext2int mapping for wind farms
% first apply mapping to all generators

order = mpc.order.gen.e2i;
% keep values corresponding to curtailable P
bindum = false(size(order));
for i=1:length(order)
    if ismember(order(i),optns.gen.curtailableP)
        bindum(i) = true;
    end
end
windorder = order(bindum);
windorder = windorder-min(windorder)+1;


if size(optns.gen.windScenarios,1) > 1
    windScenarios = optns.gen.windScenarios(windorder,:);
else
    % same production for all wind farms
    windScenarios = optns.gen.windScenarios;
end
%[windScenarios mpc.gen(mpc.gen2(:,PTYPE)==PCUR,PG)]

if isempty(optns.gen.windProbabilities)
    optns.gen.windProbabilities = 1/nScenarios*ones(1,nScenarios);
end
if optns.gen.useWindScenarios
   list = repmat(list,nScenarios,1);

   % add column with wind scenarios
   if size(list,2) < CONT_WIND_SCENARIO
       list = [list zeros(size(list,1),1)];
   end
   for i=1:nScenarios
        list(1+(i-1)*N:i*N,CONT_WIND_SCENARIO) = i;
   end
   % update probabilities
   for i=1:size(list,1)
        list(i,CONT_PROB) = list(i,CONT_PROB) * ...
                optns.gen.windProbabilities(list(i,CONT_WIND_SCENARIO));
   end
   assert(abs(sum(list(:,CONT_PROB))-1)<1e-10,'Probabilities not equal to 1');
end

list = [zeros(1,size(list,2));list]; % filler row for base case

% update number of contingencies
N = size(list,1);


nb = size(mpc.bus,1);
nl = size(mpc.branch,1);
ng = size(mpc.gen,1);

[bus,bus2,branch,gen2,gen] = deal(mpc.bus,mpc.bus2,mpc.branch,mpc.gen2,mpc.gen);


idxCurtail = find(gen2(:,PTYPE) == PCUR); % note e2i
nCurtail = length(idxCurtail);

% find lines with constraints
if optns.branch.limit
    l2 = branch(:,RATE_A) ~= 0;
    nl2 = sum(l2);
else
    nl2 = 0;
end
    

% find buses in load increas area
varload_idx = find(bus2(:,LOAD_INCREASE_AREA));
baseload = bus(:,[PD QD]);
busVarN = length(varload_idx); % number of variable buses

% containers for admittance matrices
Ybus = sparse(N*nb,nb);
Yf = sparse(N*nl,nb);
Yt = sparse(N*nl,nb);
load = sparse(N*busVarN,2); % container for load [PD QD]
wind = zeros(N*nCurtail,1); % container for wind scenarios
prob = zeros(N,1); % probabilities of scenarios
activeLines = true(nl,N);
constrainedActiveLines = false(nl,N);
activeGenerators = true(ng,N);
counter = 0;
for i=1:N
    
    inl = nl;
    impc = mpc;
    if i > 1 % contingency
        
        % apply load increase
        load_mult = list(i,CONT_LOAD_MULT);
        
        % apply wind scenario
        if list(i,CONT_WIND_SCENARIO)
%             wind(1+(i-1)*nCurtail:i*nCurtail) = optns.gen.windScenarios(...
%                  :,... % Convert to INTERNAL indexing
%                 list(i,CONT_WIND_SCENARIO) );
            wind(1+(i-1)*nCurtail:i*nCurtail) = windScenarios(:, list(i,CONT_WIND_SCENARIO) );
        else % use base case values
             wind(1+(i-1)*nCurtail:i*nCurtail) = gen(idxCurtail,PG);
        end
        
        % apply contingency
        switch list(i,CONT_TYPE) %% NOTE: mpc is internal indexing, contingency list external
            case NO_FAULT % no change
                
                
            case TRIP_LINE % implement other contingency types
                lmult = list(i,CONT_LINE_MULT);
                lidx = list(i,CONT_IDX);
                if lmult > 0
                    % increase X,R, decrease B
                    impc.branch(lidx,[BR_R BR_X BR_B]) = ...
                        [impc.branch(lidx,[BR_R BR_X])/lmult ...
                        impc.branch(lidx,BR_B)*lmult ];
                else % remove line
                    impc.branch(lidx,:) = [];
                    activeLines(lidx,i) = false;
                    inl = nl -1;
                end
                
            case TRIP_GEN
                % note: list uses external indexing
				genidx = mpc.order.gen.i2e(list(i,CONT_IDX));
				activeGenerators(genidx,i) = false;
				
			
        end
    else % base case
        load_mult = 1;
        % enter original PG into wind vector
        wind(1:nCurtail) = mpc.gen(mpc.gen2(:,PTYPE)==PCUR,PG);
    end
    
    if nl2 % prepare for line flow constraints
		
        %%- NOTE: may want to add a counter here to get the total number of
        % branch flow constraints, taking into account that tripped lines
        % in contingency cases will not be included in the constraints.
        % Otherwise dh,h will have to change size in the loop in vscopf_g
		constrainedActiveLines(:,i) = (branch(:,RATE_A) ~= 0) .* (activeLines(:,i));
    end
    
    [iYbus,iYf,iYt] = makeYbus(impc.baseMVA,impc.bus,impc.branch);
    % store contingency admittance matrix
    Ybus(1+(i-1)*nb:i*nb,:) = iYbus;
    Yf(counter+1:counter+inl,:) = iYf;
    Yt(counter+1:counter+inl,:) = iYt;
    
    
    % store load
        load(1+busVarN*(i-1):busVarN*i,:) = ...
            baseload(varload_idx,:) * load_mult;

	counter = counter + inl;		
end
nConstrainedActiveLines = sum(constrainedActiveLines,1);
nActiveLines = sum(activeLines,1);
nActiveGenerators = sum(activeGenerators,1);


c = struct(); % create struct
[c.list, c.N, c.Ybus, c.Yf, c.Yt, c.load, c.activeLines, c.nActiveLines, c.activeGenerators, c.nActiveGenerators, c.constrainedActiveLines, c.nConstrainedActiveLines, c.wind]...
	= deal(list,N,Ybus,Yf,Yt,load, activeLines, nActiveLines, activeGenerators, nActiveGenerators, constrainedActiveLines, nConstrainedActiveLines,wind); % assign elements to struct
[c.probabilities, c.nScenarios, c.nContingencies] = deal(list(:,CONT_PROB), nScenarios, nContingencies);
mpc.contingencies = c; % save struct in mpc
end
