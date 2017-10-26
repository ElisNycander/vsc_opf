function mpc = setup_contingencies(mpc,optns)

define_constants;

list = eval(optns.contingencyFile);
N = size(list,1);
nb = size(mpc.bus,1);
[bus,bus2] = deal(mpc.bus,mpc.bus2);

% find buses in load increas area
varload_idx = find(bus2(:,LOAD_INCREASE_AREA));
baseload = bus(:,[PD QD]);
busVarN = length(varload_idx); % number of variable buses

% containers for admittance matrices
Ybus = sparse((N+1)*nb,nb);
Yf = sparse((N+1)*nb,nb);
Yt = sparse((N+1)*nb,nb);
load = sparse((N+1)*busVarN,2); % container for load [PD QD]

for i=1:N+1
    
    
    impc = mpc;
    if i > 1 % 1 = base case
        sm = mpc.stabilityMargin;
        % apply contingency
        switch list(CONT_TYPE,i-1) %% NOTE: mpc is internal indexing, contingency list external
            case 0 % no change
                
                
            case 1 % implement other contingency types
        end
    else
        sm = 0;
    end
    [iYbus,iYf,iYt] = makeYbus(impc.baseMVA,impc.bus,impc.branch);
    % store contingency admittance matrix
    Ybus(1+(i-1)*nb:(i)*nb,:) = iYbus;
    Yf(1+(i-1)*nb:(i)*nb,:) = iYf;
    Yt(1+(i-1)*nb:(i)*nb,:) = iYt;
    
    % store load
        load(1+busVarN*(i-1):busVarN*i,:) = ...
            baseload(varload_idx,:) * (1 + sm);

    
end


c = struct(); % create struct
[c.list, c.N, c.Ybus, c.Yf, c.Yt, c.load] = deal(list,N,Ybus,Yf,Yt,load); % assign elements to struct
mpc.contingencies = c; % save struct in mpc
end
