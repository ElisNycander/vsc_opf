function [success, mpcArray] = verify_solutions(results,om,optns)
% Verify solutions by running power flows. 
% The values from the tables in results are entered into the mpc struct for
% each case, and the power flow should in the first iteration, except when
% it needs to make conversions from PV buses to PQ buses
% Note that it would be possible to determine from the lagrange multipliers
% which buses should be converted as well.

%clear;

%load('verify_solutions.mat'); 



define_constants;

mpc = get_mpc(om); % Note: mpc here uses internal indexing
[cs,bus2,baseMVA,version,gen2] = deal(mpc.contingencies,mpc.bus2,mpc.baseMVA,mpc.version,mpc.gen2);
[mpopt] = deal(optns.mpopt);
[N,list,nActiveLines, activeLines] = deal(cs.N,cs.list,cs.nActiveLines,cs.activeLines);
idxLoadIncreaseBuses = bus2(:,LOAD_INCREASE_AREA) == 1;
nLoadIncreaseBuses = sum(idxLoadIncreaseBuses);
 
idxPQFactor = gen2(:,PQ_FACTOR) ~= 0;
PQFactors = gen2(idxPQFactor,PQ_FACTOR);

% do pfs quietly
mpopt.out.all = 0;
mpopt.verbose = 0;

ycounter = 0;
success = ones(1,N);
for i=1:N
   
   
    [bus,gen,branch] = deal(mpc.bus,mpc.gen,mpc.branch);
 
	%% bus matrix
   % enter voltages into bus matrix
    Vm = table2array(results.Vm(:,i+3));
	Va = table2array(results.Va(:,i+1));
    bus(:,VM) = Vm;
    bus(:,VA) = Va;
	
	% load
	bus(idxLoadIncreaseBuses,[PD QD]) = cs.load(1+(i-1)*nLoadIncreaseBuses:i*nLoadIncreaseBuses,:);
    
   %% gen matrix	
    % deactivate generators

    %gen = gen(cs.activeGenerators(:,i),:); % remove inactive generators
	
	Pg = table2array(results.Pg(:,i+4)); % Note: Pg and Qg uses external indexing
	%Pg = Pg(~isnan(Pg)); % remove NaN for inactive generators
    % convert to internal indexing
    Pg = Pg(mpc.order.gen.e2i,:);
    
    
	Qg = table2array(results.Qg(:,i+4));
	%Qg = Qg(~isnan(Pg));
    % convert to internal indexing
    Qg = Qg(mpc.order.gen.e2i,:);
    
	gen(:,[PG QG]) = [Pg Qg];
    
    %% modified Q-limits
    % Note: Those generators with PQ-constraints effectively have reduced Q
    % limits. These must be carried over to the new mpc struct to get the
    % same results. Note that PQ-constraints are only applied to
    % contingency cases
    if (i>1) && optns.gen.usePQConstraints
        gen(idxPQFactor,[QMIN QMAX]) = [-PQFactors.*abs(Pg(idxPQFactor)) ...
            PQFactors.*abs(Pg(idxPQFactor)) ];
    end
        
    % remove inactive generators
    idxInactive = isnan(gen(:,PG));
    gen(idxInactive,:) = [];
    
    % enter bus voltages into gen
    for ii=1:size(gen,1)
        gen(ii,VG) = bus(gen(ii,GEN_BUS)==bus(:,BUS_I),VM);
    end
    
	
	%% branch matrix
	% change impedance
	if list(i,CONT_TYPE) == 1 && list(i,CONT_LINE_MULT) ~= 0
		lineNr = list(i,CONT_IDX);
		branch(lineNr,[BR_R BR_X]) = branch(lineNr,[BR_R BR_X]) / list(i,CONT_LINE_MULT);
		branch(lineNr,BR_B) = branch(lineNr,BR_B) * list(i,CONT_LINE_MULT);
	end
	% remove tripped lines
	branch = branch(cs.activeLines(:,i),:);
	
    
	impc = struct();
    %mpopt.pf.enforce_q_lims = 0;
	[impc.bus,impc.gen,impc.branch,impc.baseMVA,impc.version] = deal(bus,gen,branch,baseMVA,version);
	impc = runpf(impc,mpopt);
	if impc.iterations > 0
        V0 = Vm.*exp(1j*pi/180*Va);
        V1 = impc.bus(:,VM).*exp(1j*pi/180*impc.bus(:,VA));
        diff = sum(abs(V0-V1));
        if diff > mpopt.pf.tol
            success(i) = 0;
        end
    end
    
    mpcArray(i) = impc;
	
	ycounter = ycounter + nActiveLines(i);
end

