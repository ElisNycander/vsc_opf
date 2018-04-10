function cs = maximum_transfer_vs_qlim(mpc,sbus,sQlim,caseStr)
% Find transfer/loadability limit as function of reactive power capability
% in given bus

%% for script
% clear;
% close all;
% 
% % load case data
% mpc = Nordic32;
% mpc = define_areas(mpc);
% 
% caseStr = 'Nordic 32 unmodified';
% 
% sbus = 4032;
% %sQlim = 0:50:1e3; % range of Q-values
% sQlim = 350;


%% function
%display('%%%%%%%%%%%%%%%%%%%% NEW RUN %%%%%%%%%%%%%%%%%%%%%%%%');
desc = caseStr;
findBaseTransfer = 1;


define_constants;
%define_constants_N32;
spica_settings;

plotCurve = 0;

detectNoseBus = [4044 4045 4046]; % buses used for detecting if we move to an unstable manifold
generation_increase = 1; % Also increase generation, if 0 only slack bus compensates increase
target_transfer_increase = 100;

ng = size(mpc.gen,1);


plot_buses = [  
    4021
%               4071 
%               4072
               4011
%               4012
%               4022
%               4021
               4031
               4032
              4041
              4042
              4043
              4044
              4045
              4046
              4047
              4051
%               4061
%               4062
%               4063
%                 41
%                 42
%                 43
%                 46
%                 47
%                 51
%                 61
%                 62
%                 63
              ];

%% options

mpopt  =  mpoption('out.all',  1,  'verbose',  2);

mpopt = mpoption(mpopt, 'out.sys_sum', 1,  ...
                        'out.area_sum', 1, ...
                        'out.bus', 1, ...
                        'out.branch', 0, ...
                        'out.gen', 1, ...
                        'out.force', 0 );

% pf options
mpopt = mpoption(mpopt, 'pf.nr.max_it', 20);
mpopt = mpoption(mpopt, 'pf.enforce_q_lims', 1);

% cpf options
mpopt = mpoption(mpopt, 'cpf.stop_at', 'NOSE', 'cpf.step', 0.1); 
mpopt = mpoption(mpopt, 'cpf.plot.level', 0);
mpopt = mpoption(mpopt, 'cpf.plot.bus', plot_buses);
mpopt = mpoption(mpopt, 'cpf.parameterization', 2);
mpopt = mpoption(mpopt, 'cpf.enforce_q_lims',1);
mpopt = mpoption(mpopt, 'cpf.adapt_step',0);
mpopt = mpoption(mpopt, 'cpf.user_callback','n32_callback');


%% add generator with p=PG to bus
sbus_idx = find( mpc.bus(:,BUS_I) == sbus );

% set bus type to PV
mpc.bus( sbus_idx,BUS_TYPE ) = 2;

% find other generators at bus
sGidx = find( mpc.gen( :,GEN_BUS ) == sbus);
% current Q-capability at bus
sQMAX0 = sum( mpc.gen( sGidx,QMAX ) );
sQMIN0 = sum( mpc.gen( sGidx,QMIN ) );
sPG0 = sum( mpc.gen( sGidx,PG ) );

% use voltage of old generators at bus
Vg = min( mpc.gen( sGidx,VG ) );
if isempty(Vg)
    Vg = 1;
end
    
% remove other generators
mpc.gen(sGidx,:) = [];
% add generation as negative load
%mpc.bus( sbus_idx,PD ) = mpc.bus( sbus_idx,PD ) - sPG0;


newgen = [
    sbus   sPG0 0   0      0       Vg   100     1       1e3     0
];
mpc.gen = [mpc.gen;
    newgen zeros( 1,size(mpc.gen,2)-size( newgen,2 ) )
];

sgen = size( mpc.gen,1 );

%% lower Qlims

mpc.gen( :,QMIN ) = -1e3;

%% function
plot_flows = 0;


% clear file
fileID = fopen('nosepf.txt','w');
fclose(fileID);
fileID = fopen('cpf.txt','w');
fclose(fileID);


%% indices of buses to plot
pi = false(1,size(mpc.bus,1));
for i=1:length(plot_buses)
    idx = bus_idx(plot_buses(i),mpc);
    pi(idx) = true;
end

% generators with generation increase
idxGenIncrease = [];
for i=1:size(spica_gen_increase,1)
    bus = spica_gen_increase(i,1)
    find(mpc.gen(:,GEN_BUS)==bus)
    idxGenIncrease = [idxGenIncrease find(mpc.gen(:,GEN_BUS)==bus)]
end

% buses with load increase
idxLoadIncrease = [];
for i=1:size(spica_load_increase,1)
    bus = spica_load_increase(i,1);
    idxLoadIncrease = [idxLoadIncrease find(mpc.bus(:,1)==bus)];
end

% find branches in bottleneck
nSpicaBranch = size(spica_bottleneck,1);
idxSpicaBranch = [];
for ii=1:nSpicaBranch
    bus1 = spica_bottleneck(ii,1);
    bus2 = spica_bottleneck(ii,2);
    
    idxSpicaBranch = [idxSpicaBranch line_idx(bus1,bus2,mpc)];
end

N = length(sQlim);

cpf_strs = cell(1,N);
cpf_flag = zeros(1,N);
cpf_max_transfer = zeros(1,N);

for i=1:N
    %% base case
    mpcb = mpc;
    
    % change Q capability
    mpcb.gen(sgen,QMAX) = sQlim(i);
    
    
    %% solve base case
    try mpcb = runpf(mpcb,mpopt);
    catch
        error(['Base case for scenario nr %d not solvable: ' mpcb.desc],k);
    end
    
    %% target case
    mpct = mpcb;
    
    PGtot = sum(mpcb.gen(idxGenIncrease,PG));
    % increase generation
    if generation_increase
        mpct.gen(idxGenIncrease,PG) = mpcb.gen(idxGenIncrease,PG) + mpcb.gen(idxGenIncrease,PG)*target_transfer_increase/PGtot;
    end
    
    PDtot = sum(mpcb.bus(idxLoadIncrease,PD));
    % increase load
    mpct.bus(idxLoadIncrease,[PD QD]) = mpcb.bus(idxLoadIncrease,[PD QD]) + mpcb.bus(idxLoadIncrease,[PD QD])*target_transfer_increase/PDtot;
    
    
    %% run cpf
    [mpcc, flag] = runcpf(mpcb,mpct,mpopt,'cpf.txt');
    
    
    %% Check if voltages are increasing in the end - have transitioned 
    %% to unstable manifold
    % exclude current bus from check (only PQ-buses)
    noseBuses = setdiff(detectNoseBus,sbus);
    dbus_idx = [];
    for ii=1:length(noseBuses)
        dbus_idx = [dbus_idx find( mpc.bus( :,BUS_I) == noseBuses(ii) )];   
    end
    
    V = mpcc.cpf.V;
    Vm = abs( V(dbus_idx,:) );
    
    [~,I] = min( Vm,[],2);
    
    flag = 1;
    for ii=1:length(I)-1
        if I(ii) ~= I(ii+1)
            flag = 0;
        end
    end
    if ~flag
        warning('Minima for voltage correspond to different lambda at different buses');
        idx = mode(I,'minmode');
    else
        idx = I(1);
     end

    Vnose = V(:,idx);
    
    % calculate transfers for this case
    mpc_int = ext2int(mpc);
    [~,YF,YT] = makeYbus( mpc_int );
    
    iSf = Vnose(mpc_int.branch(:, F_BUS)) .* conj(YF * Vnose);  %% complex power injected at "from" bus (p.u.)
    iSt = Vnose(mpc_int.branch(:, T_BUS)) .* conj(YT * Vnose);  %% complex power injected at "to" bus (p.u.)
    
    % put transfers into mpcc
    mpcc.branch( mpcc.order.branch.status.on, [PF PT QF QT] ) = mpc.baseMVA * [
        real(iSf) real(iSt) imag(iSf) imag(iSt) 
    ];    
   %% find total transmission through bottleneck
    
    % sum transfer over bottleneck
    transfers = zeros(3,4);
    for ii=1:length(idxSpicaBranch)
        if idxSpicaBranch > 0
            transfers( 1,: ) = transfers( 1,: ) ...
                + mpcb.branch( idxSpicaBranch(ii), [PF PT QF QT] );
            transfers( 2,: ) = transfers( 2,: ) ...
                + mpct.branch( idxSpicaBranch(ii), [PF PT QF QT] );
            transfers( 3,: ) = transfers( 3,: ) ...
                + mpcc.branch( idxSpicaBranch(ii), [PF PT QF QT] );
        else
            transfers( 1,: ) = transfers( 1,: ) ...
                + mpcb.branch( idxSpicaBranch(ii), [PT PF QT QF] );
            transfers( 2,: ) = transfers( 2,: ) ...
                + mpct.branch( idxSpicaBranch(ii), [PT PF QT QF] );
            transfers( 3,: ) = transfers( 3,: ) ...
                + mpcc.branch( idxSpicaBranch(ii), [PT PF QT QF] );
        end
        % Note: Target case not solved, gives same values as base case
        
    end
    
    mpcb.transfer = transfers(1,:);
    mpcc.transfer = transfers(3,:);
    
    
    %% Create string with information on node conversion
    
    cpf_strs{i} = create_cpf_str(desc,mpcc,mpcb,mpc);
    cpf_flag(i) = flag;
    cpf_max_transfer(i) = mpcc.transfer(1);
    %disp(i);
end

%% try to find sQMAX0

if findBaseTransfer
    sidx = find( sQlim == sQMAX0 );
    if ~isempty(sidx)
        tx0 = cpf_max_transfer(sidx);
    else % interpolate between closest values
        sidx = find( sQlim > sQMAX0,1 );
        
        tx0 = cpf_max_transfer(sidx-1) + (cpf_max_transfer(sidx) - cpf_max_transfer(sidx-1)) ...
            / ( sQlim(sidx) - sQlim(sidx-1) ) * ( sQMAX0 - sQlim(sidx-1) );
    end
end

if plotCurve
    figure;
    PH1 = plot(sQlim, cpf_max_transfer);
    PH1.Marker = 'o';
    grid on;
    hold on;
    title(['Transfer capacity vs QLIM at bus ' ...
        num2str(sbus) sprintf('\n') desc]);
    xlabel('Qmax (MW)');
    ylabel('Transfer (MW)');
    
    LH1 = line([sQMAX0 sQMAX0],ylim);
    LH1.Color = [0 0 0];
    LH1.LineStyle = '--';
    if findBaseTransfer
        plot( sQMAX0, tx0 ,'*');
    end
end

cs = struct();

if findBaseTransfer
[cs.maxTransfer, cs.str, cs.baseTransfer, cs.baseQ] = ...
    deal(cpf_max_transfer,cpf_strs,tx0,sQMAX0);
else
    [cs.maxTransfer, cs.str, cs.baseQ] = ...
    deal(cpf_max_transfer,cpf_strs,sQMAX0);
end