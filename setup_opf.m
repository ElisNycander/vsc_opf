function om = setup_opf(mpc,optns)
%OPF  Constructs an OPF model object from a MATPOWER case struct.
%   OM = OPF_SETUP(MPC, MPOPT)
%

define_constants;


%% assign variables for simplicity
[mpopt,gen2,cs] = deal(optns.mpopt,mpc.gen2,mpc.contingencies);

%% data dimensions
nb   = size(mpc.bus, 1);    %% number of buses
nl   = size(mpc.branch, 1); %% number of branches
ng   = size(mpc.gen, 1);    %% number of dispatchable injections
nc   = cs.N; %% number of contingencies
nPfix = sum(mpc.gen2(:,PFIX)); %% number of fixed generators
nQfix = sum(mpc.gen2(:,QFIX));
idxPfix = mpc.gen2(:,PFIX) == 1;
idxQfix = mpc.gen2(:,QFIX) == 1;

if optns.branch.limit
    idxConstrainedLines = mpc.branch(:,RATE_A) ~= 0;
else
    idxConstrainedLines = zeros(size(mpc.branch,1),1);
end
nConstrainedLines = sum(idxConstrainedLines);

%% create (read-only) copies of individual fields for convenience
[baseMVA, bus, gen, branch, gencost, Au, lbu, ubu, mpopt, ...
    N, fparm, H, Cw, z0, zl, zu, userfcn] = opf_args(mpc, mpopt);



%% set up initial variables and bounds
Va   = bus(:, VA) * (pi/180);
Vm   = bus(:, VM);
Pg   = gen(:, PG) / baseMVA;
Qg   = gen(:, QG) / baseMVA;
Pmin = gen(:, PMIN) / baseMVA;
Pmax = gen(:, PMAX) / baseMVA;
Qmin = gen(:, QMIN) / baseMVA;
Qmax = gen(:, QMAX) / baseMVA;

% initial variables for contingencies
% Vac = zeros(size(Va));
% Vmc = ones(size(Vm));
% Qgc = zeros(size(Qg));
Vac = Va; Vmc = Vm; Qgc = Qg;

%% contingency variables
ngc_idx = ~logical(gen2(:,PFIX));
% only some generator may vary active power
Pg_c = Pg(ngc_idx);
Pmin_c = Pmin(ngc_idx);
Pmax_c = Pmax(ngc_idx);

%% warn if there is more than one reference bus
refs = find(bus(:, BUS_TYPE) == REF);
if length(refs) > 1 && mpopt.verbose > 0
  errstr = ['\nopf_setup: Warning: Multiple reference buses.\n', ...
              '           For a system with islands, a reference bus in each island\n', ...
              '           may help convergence, but in a fully connected system such\n', ...
              '           a situation is probably not reasonable.\n\n' ];
  fprintf(errstr);
end

%% voltage angle reference constraints
Vau = Inf(nb, 1);
Val = -Vau;
Vau(refs) = Va(refs);
Val(refs) = Va(refs);


%% construct OPF model object
om = opf_model(mpc);

  % contingency variables
  % note that for each contingency we have a new set of
  % Vm,Va,Pg(ngc_idx),Qg
  for i=1:nc
      
      if i == 1 % base case
          % base case
          om = add_vars(om, 'Va', nb, Va, Val, Vau);
          om = add_vars(om, 'Vm', nb, Vm, bus(:, VMIN), bus(:, VMAX));
          om = add_vars(om, 'Pg', ng, Pg, Pmin, Pmax);
          om = add_vars(om, 'Qg', ng, Qg, Qmin, Qmax);
          
          om = add_constraints(om, 'Pmis', nb, 'nonlinear');
          om = add_constraints(om, 'Qmis', nb, 'nonlinear');
          
          om = add_constraints(om, 'Sf', nConstrainedLines, 'nonlinear');
          om = add_constraints(om, 'St', nConstrainedLines, 'nonlinear');
          
      else
          
          stringIdx = num2str(i);
          
          om = add_vars(om, ['Va' stringIdx], nb, ...
              repmat(Vac,1,1), ...
              repmat(Val,1,1), ...
              repmat(Vau,1,1) ...
              );
          om = add_vars(om, ['Vm' stringIdx], nb, ...
              repmat(Vmc,1,1), ...
              repmat(bus(:,VMIN),1,1), ...
              repmat(bus(:,VMAX),1,1) ...
              );
          % note: number of generators may change if contingency includes
          % outage of generator
          
          idxQvar = and(cs.activeGenerators(:,i), ~idxQfix);
          idxPvar = and(cs.activeGenerators(:,i), ~idxPfix);
          nQvar = sum(idxQvar);
          nPvar = sum(idxPvar);

          iPg = Pg(idxPvar);
          iPmin = Pmin(idxPvar);
          iPmax = Pmax(idxPvar);
          iQg = Qgc(idxQvar);
          iQmin = Qmin(idxQvar);
          iQmax = Qmax(idxQvar);

          %ngc  = ng - nfix; %% number of generators that may increase generation in stressed cases
          om = add_vars(om, ['Pg' stringIdx], nPvar, iPg, iPmin, iPmax);
          om = add_vars(om, ['Qg' stringIdx], nQvar, iQg, iQmin, iQmax);
          
          %%- DONT ADD CONSTRAINTS AS THESE ARE NOT USED
          % power balance
          om = add_constraints(om, ['Pmis' stringIdx], nb, 'nonlinear');
          om = add_constraints(om, ['Qmis' stringIdx], nb, 'nonlinear');
          % branch flows - may be less than nConstrainedLines for line outages
          %nConstrainedActiveLines = sum(idxConstrainedLines .* cs.activeLines(:,i));
          om = add_constraints(om, ['Sf' stringIdx], cs.nConstrainedActiveLines(i), 'nonlinear');
          om = add_constraints(om, ['St' stringIdx], cs.nConstrainedActiveLines(i), 'nonlinear');
      end
  end
