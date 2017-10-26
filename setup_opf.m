function om = setup_opf(mpc,optns)
%OPF  Constructs an OPF model object from a MATPOWER case struct.
%   OM = OPF_SETUP(MPC, MPOPT)
%

define_constants;


%% assign variables for simplicity
[mpopt,gen2,contingencies] = deal(optns.mpopt,mpc.gen2,mpc.contingencies);

%% data dimensions
nb   = size(mpc.bus, 1);    %% number of buses
nl   = size(mpc.branch, 1); %% number of branches
ng   = size(mpc.gen, 1);    %% number of dispatchable injections
nc   = contingencies.N; %% number of contingencies
nfix = sum(mpc.gen2(:,PFIX)); %% number of fixed generators
ngc  = ng - nfix; %% number of generators that may increase generation in stressed cases


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

% base case
om = add_vars(om, 'Va', nb, Va, Val, Vau);
om = add_vars(om, 'Vm', nb, Vm, bus(:, VMIN), bus(:, VMAX));
om = add_vars(om, 'Pg', ng, Pg, Pmin, Pmax);
om = add_vars(om, 'Qg', ng, Qg, Qmin, Qmax);

om = add_constraints(om, 'Pmis', nb, 'nonlinear');
om = add_constraints(om, 'Qmis', nb, 'nonlinear');

om = add_constraints(om, 'Sf', nl, 'nonlinear');
om = add_constraints(om, 'St', nl, 'nonlinear');

  % contingency variables
  % note that for each contingency we have a new set of
  % Vm,Va,Pg(ngc_idx),Qg
  for i=1:nc
      om = add_vars(om, ['Va' num2str(i)], nb, ...
          repmat(Va,nc,1), ...
          repmat(Val,nc,1), ...
          repmat(Vau,nc,1) ...
          );
      om = add_vars(om, ['Vm' num2str(i)], nb, ...
          repmat(Vm,nc,1), ...
          repmat(bus(:,VMIN),nc,1), ...
          repmat(bus(:,VMAX),nc,1) ...
          );
      % note: number of generators may change if contingency includes
      % outage of generator
      om = add_vars(om, ['Pg' num2str(i)], ngc, ...
          repmat(Pg_c,nc,1), ...
          repmat(Pmin_c,nc,1), ...
          repmat(Pmax_c,nc,1) ...
          );
      om = add_vars(om, ['Qg' num2str(i)], ng, ...
          repmat(Qg,nc,1), ...
          repmat(Qmin,nc,1), ...
          repmat(Qmax,nc,1) ...
          );
      

      % power balance
      om = add_constraints(om, ['Pmis' num2str(i)], nb, 'nonlinear');
      om = add_constraints(om, ['Qmis' num2str(i)], nb, 'nonlinear');
      % branch flows
      om = add_constraints(om, ['Sf' num2str(i)], nl, 'nonlinear');
      om = add_constraints(om, ['St' num2str(i)], nl, 'nonlinear');
  end
