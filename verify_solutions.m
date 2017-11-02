function [success, mpcArray] = verify_solutions(results,om,optns)

%clear;

%load('verify_solutions.mat'); 



define_constants;

mpc = get_mpc(om);
[cs,bus2,baseMVA,version] = deal(mpc.contingencies,mpc.bus2,mpc.baseMVA,mpc.version);
[mpopt] = deal(optns.mpopt);
[N,list,nActiveLines, activeLines] = deal(cs.N,cs.list,cs.nActiveLines,cs.activeLines);
idxLoadIncreaseBuses = bus2(:,LOAD_INCREASE_AREA) == 1;
nLoadIncreaseBuses = sum(idxLoadIncreaseBuses);
 
% do pfs quietly
mpopt.out.all = 0;
mpopt.verbose = 0;

ycounter = 0;
success = 1;
for i=1:N
   
   
    [bus,gen,branch] = deal(mpc.bus,mpc.gen,mpc.branch);
 
	%% bus matrix
   % enter voltages into bus matrix
    Vm = table2array(results.Vm(:,i+1));
	Va = table2array(results.Va(:,i+1));
    bus(:,VM) = Vm;
    bus(:,VA) = Va;
	
	% load
	bus(idxLoadIncreaseBuses,[PD QD]) = cs.load(1+(i-1)*nLoadIncreaseBuses:i*nLoadIncreaseBuses,:);
    
   %% gen matrix	
    % deactivate generators
    gen = gen(cs.activeGenerators(:,i),:);
	
	Pg = table2array(results.Pg(:,i+1));
	Pg = Pg(~isnan(Pg));
	Qg = table2array(results.Qg(:,i+1));
	Qg = Qg(~isnan(Pg));
	gen(:,[PG QG]) = [Pg Qg];
	
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
            success = 0;
        end
    end
    
    mpcArray(i) = impc;
	
	ycounter = ycounter + nActiveLines(i);
end

