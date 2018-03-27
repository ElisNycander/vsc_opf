function [success, mismatchMatrix] = verify_solutions_sbus(restab,om,optns)
% Verify solutions by running power flows. 
% The values from the tables in restab are entered into the mpc struct for
% each case, and the power flow should in the first iteration, except when
% it needs to make conversions from PV buses to PQ buses
% Note that it would be possible to determine from the lagrange multipliers
% which buses should be converted as well.

%% for script
% clear;
% 
% load('verify_solutions.mat'); 
% load('N32_ybus.mat');
% optns.verifyTol = 1e-8;

%% function
define_constants;

mpc = get_mpc(om); % Note: mpc here uses internal indexing

[cs,bus2,baseMVA,version,gen2] = deal(mpc.contingencies,mpc.bus2,mpc.baseMVA,mpc.version,mpc.gen2);
%[mpopt] = deal(optns.mpopt);
[N,list,nActiveLines, activeLines] = deal(cs.N,cs.list,cs.nActiveLines,cs.activeLines);
idxLoadIncreaseBuses = bus2(:,LOAD_INCREASE_AREA) == 1;
nLoadIncreaseBuses = sum(idxLoadIncreaseBuses);
 

ycounter = 0;
success = ones(1,N);
mismatchMatrix = zeros( 2 * size( mpc.bus,1 ),N );
for i=1:N
   
   % Note: bus, gen - INTERNAL ordering
   
    [bus,gen,branch] = deal(mpc.bus,mpc.gen,mpc.branch);
    if i == 1
        sidx = '';
    else
        sidx = num2str(i);
    end
    
	%% bus matrix
   % enter voltages into bus matrix
    Vm = table2array( restab.Vm( :,['VM' sidx] ) );
	Va = table2array( restab.Va( :,['VA' sidx] ) );
    bus(:,VM) = Vm;
    bus(:,VA) = Va;
	% complex voltages
    V = Vm .* exp( 1j * Va/180*pi );
    
	% load
	bus(idxLoadIncreaseBuses,[PD QD]) = cs.load(1+(i-1)*nLoadIncreaseBuses:i*nLoadIncreaseBuses,:);
    
   %% gen matrix	
   
	Pg = table2array(restab.Pg(:,['PG' sidx])); % Note: Pg and Qg uses external indexing
    % convert to internal indexing
    Pg = Pg(mpc.order.gen.e2i,:);
    Pg( isnan( Pg ) ) = 0;
    
	Qg = table2array(restab.Qg(:,['QG' sidx]));
    % convert to internal indexing
    Qg = Qg(mpc.order.gen.e2i,:);
    Qg( isnan( Qg ) ) = 0;
    
	gen(:,[PG QG]) = [Pg Qg];
    	
	%% branch matrix
	% change impedance
	if list(i,CONT_TYPE) == 1 && list(i,CONT_LINE_MULT) ~= 0
		lineNr = list(i,CONT_IDX);
		branch(lineNr,[BR_R BR_X]) = branch(lineNr,[BR_R BR_X]) / list(i,CONT_LINE_MULT);
		branch(lineNr,BR_B) = branch(lineNr,BR_B) * list(i,CONT_LINE_MULT);
	end
	% remove tripped lines
	branch = branch(cs.activeLines(:,i),:);
	
    Sbus = makeSbus( baseMVA, bus, gen );
    Ybus = makeYbus( baseMVA, bus, branch );
    mis = V .* conj(Ybus * V) - Sbus;
    
    mismatch = abs( [ real(mis) ; imag(mis) ] );
    mismatchMatrix( :,i ) = mismatch;
    
    if ~ isempty( find( mismatch > optns.verifyTol,1 ) )
        success(i) = 0;
    end
	
	ycounter = ycounter + nActiveLines(i);
end

