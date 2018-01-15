function mpc = setup_contingencies(mpc,optns)

define_constants;

list = eval(optns.contingencyFile);
list = [zeros(1,size(list,2));list]; % filler row for base case
if isempty(list) % minimal list if there are no contigencies
    list = zeros(1,4);  
end
N = size(list,1);

nb = size(mpc.bus,1);
nl = size(mpc.branch,1);
ng = size(mpc.gen,1);
[bus,bus2,branch,gen2] = deal(mpc.bus,mpc.bus2,mpc.branch,mpc.gen2);

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
activeLines = true(nl,N);
constrainedActiveLines = false(nl,N);
activeGenerators = true(ng,N);
counter = 0;
for i=1:N
    
    inl = nl;
    impc = mpc;
    if i > 1 % contingency
        load_mult = list(i,CONT_LOAD_MULT);
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
				genidx = list(i,CONT_IDX);
				activeGenerators(genidx,i) = false;
				
			
        end
    else % base case
        load_mult = 1;
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
[c.list, c.N, c.Ybus, c.Yf, c.Yt, c.load, c.activeLines, c.nActiveLines, c.activeGenerators, c.nActiveGenerators, c.constrainedActiveLines, c.nConstrainedActiveLines]...
	= deal(list,N,Ybus,Yf,Yt,load, activeLines, nActiveLines, activeGenerators, nActiveGenerators, constrainedActiveLines, nConstrainedActiveLines); % assign elements to struct
mpc.contingencies = c; % save struct in mpc
end
