function mpc = randomize_initial_guess(mpc,optns)
% Return solved mpc case with PG or VG assigned from rand() distribution.
% If a particular distribution results in an unsolvable power flow
% a new distribution is generated until one that is solvable is found.
%

% for script
% clear;
% load('randomize_initial_data.mat');
% optns.windScenario = 'E1';
% optns.randomizeVariable = 'V'; % P or V

%% fucntion
define_constants;
vscopf_define_constants;

% power to be allocated: load - wind power - external power

% distribute power randomly weighted by capacity (to buses excluding slack
% bus)

%Pload = sum(mpc.bus(:,PD));
%Pexternal = sum(mpc.bus( mpc.bus(:,BUS_AREA) == 4,PD));
%Pwind = sum(mpc.gen(mpc.gen2(:,PTYPE) == PCUR,PG));
%Palloc = Pload - Pwind - Pexternal;

Palloc = sum( mpc.gen( optns.gen.genArea == 1, PG ) );

% PF options
op = mpoption();
op.pf.enforce_q_lims = 1;
op.verbose = 0;
op.out.all = 0;


while 1 % try random solutions until a valid solution is found
    
    switch optns.windScenario
        case 'E1' % randomly allocate
            
            switch optns.randomizeVariable
                
                case 'P'
                    %sum( mpc.gen( optns.gen.genArea == 1, PG ) )
                    %mpc.gen( optns.gen.genArea == 1, PG )
                    mpc.gen( optns.gen.genArea == 1, PG ) = rand( sum( optns.gen.genArea == 1 ) ,1 ) .* ...
                        mpc.gen( optns.gen.genArea == 1, PMAX );
                    %mpc.gen( optns.gen.genArea == 1, PG )
                    
                    % put slack to 0
                    mpc.gen( mpc.gen(:,GEN_BUS) == mpc.bus( mpc.bus(:,BUS_TYPE) == 3,BUS_I ) ,PG) = 0;
                    
                    % scale generation to Palloc
                    Pg = sum( mpc.gen( optns.gen.genArea == 1, PG ) );
                    
                    mpc.gen( optns.gen.genArea == 1, PG ) = Palloc / Pg * mpc.gen( optns.gen.genArea == 1, PG );
                    
                    %mpc.gen( optns.gen.genArea == 1, PG )
                    %sum( mpc.gen( optns.gen.genArea == 1, PG ) )
                    
                case 'V' % randomize VG for all generators
                    
             
                    % randomize VM for all buses
                    mpc.bus( :,VM ) = mpc.bus( :,VMIN ) + ( 1/2 + rand( size(mpc.bus,1),1 ) ) .* ...
                                        ( mpc.bus( :,VMAX ) - mpc.bus( :,VMIN ) ) / 2;
                    %mpc.bus( :,VM )                
                                    
                    % put bus VM into gen VG
                    for i = 1:size( mpc.gen,1 )
                        % find bus idx
                        busIdx = find( mpc.bus( :,BUS_I ) == mpc.gen( i,GEN_BUS ) );
                        mpc.gen( i,VG ) = mpc.bus( busIdx,VM );
                    end
                        
            end
    end
    % run pf
    try
        mpco = runpf(mpc,op);
    catch
        continue;
    end
    
    if mpco.success == 1
        mpc = mpco;
        break;
    end
    
end
