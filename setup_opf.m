function om = setup_opf(mpc,optns)
%OPF  Constructs an OPF model object from a MATPOWER case struct.
%   OM = OPF_SETUP(MPC, MPOPT)
%

define_constants;


%% assign variables for simplicity
[mpopt,cs] = deal(optns.mpopt,mpc.contingencies);

%% data dimensions
nb   = size(mpc.bus, 1);    %% number of buses
nl   = size(mpc.branch, 1); %% number of branches
ng   = size(mpc.gen, 1);    %% number of dispatchable injections
nc   = cs.N; %% number of contingencies


% nFix = sum(mpc.gen2(:,PTYPE)==PFIX);
% nVar = sum(mpc.gen2(:,PTYPE)==PVAR);
% nCur = sum(mpc.gen2(:,PTYPE)==PCUR);

idxPVar = mpc.gen2(:,PTYPE)==PVAR;
idxPFix = mpc.gen2(:,PTYPE)==PFIX;
% curtailable generators
idxCurtail = mpc.gen2(:,PTYPE)==PCUR;
idxPQFactor = mpc.gen2(:,PQ_FACTOR) ~= 0;
nCurtail = sum(idxCurtail);

idxQVar = ones(ng,1);

if optns.branch.limit
    idxConstrainedLines = mpc.branch(:,RATE_A) ~= 0;
else
    idxConstrainedLines = zeros(size(mpc.branch,1),1);
end
nConstrainedLines = sum(idxConstrainedLines);

%% create (read-only) copies of individual fields for convenience
% [baseMVA, bus, gen, branch, gencost, Au, lbu, ubu, mpopt, ...
%     N, fparm, H, Cw, z0, zl, zu, userfcn] = opf_args(mpc, mpopt);
[baseMVA, bus, gen, gen2] = deal(mpc.baseMVA,mpc.bus,mpc.gen, mpc.gen2);


%% set up initial variables and bounds
Va   = bus(:, VA) * (pi/180);
Vm   = bus(:, VM);
Pg   = gen(:, PG) / baseMVA;
Qg   = gen(:, QG) / baseMVA;
Pmin = gen(:, PMIN) / baseMVA;
Pmax = gen(:, PMAX) / baseMVA;
Qmin = gen(:, QMIN) / baseMVA;
Qmax = gen(:, QMAX) / baseMVA;

idxcq = and(idxPFix,idxPQFactor);
%change Q-lmits for fixed generators, with non-zero PQ-factor
if ~optns.gen.optimizeBaseP && optns.gen.usePQConstraints
    Qmin(idxcq) = -gen2(idxcq,PQ_FACTOR).*abs(Pg(idxcq));
    Qmax(idxcq) = gen2(idxcq,PQ_FACTOR).*abs(Pg(idxcq));
end

% initial variables for contingencies
% Vac = zeros(size(Va));
% Vmc = ones(size(Vm));
% Qgc = zeros(size(Qg));
Vac = Va; Vmc = Vm; Qgc = Qg;


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
          if optns.gen.optimizeBaseP % include active power for base case as optimization variables
            om = add_vars(om, 'Pg', ng, Pg, Pmin, Pmax);
          end
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
          idxActive = cs.activeGenerators(:,i);
          
          idxQvar = and(idxActive, idxQVar);
          idxPvar = and(idxActive, idxPVar);
          nQvar = sum(idxQvar);
          nPvar = sum(idxPvar);
          

          iPg = Pg(idxPvar);
          iPmin = Pmin(idxPvar);
          iPmax = Pmax(idxPvar);
          iQg = Qgc(idxQvar);
          iQmin = Qmin(idxQvar);
          iQmax = Qmax(idxQvar);
          
          % add P and Q for variable generators
          om = add_vars(om, ['Pg' stringIdx], nPvar, iPg, iPmin, iPmax);
          om = add_vars(om, ['Qg' stringIdx], nQvar, iQg, iQmin, iQmax);
          
          % add curtailment fractions
          % NOT IMPLEMENTED
          zv = zeros(nCurtail,1);
          ov = ones(nCurtail,1);
          om = add_vars(om, ['Beta' stringIdx],nCurtail,zv,zv,ov);
          
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
  
  %% linear constraints for PQ capabilities
  if optns.gen.usePQConstraints
      vv = get_idx(om);
      x0 = getv(om);
      n = length(x0);
      
      
      for i=2:nc % only for contingencies
          
          % prepare constraint matrices
          A = zeros(1,n);
          vl = zeros(1); vu = inf(1);
          
          Aidx = 1;
          
          idxActive = cs.activeGenerators(:,i);
          stringIdx = num2str(i);
          %% get variable indices
          nVarP = vv.N.(['Pg' stringIdx]);
          iPgVar = vv.i1.(['Pg' stringIdx]):vv.iN.(['Pg' stringIdx]);
          idxPvar = find(and(idxActive,idxPVar));
          assert(nVarP == length(idxPvar),'nVarp ~= length(idxPvar)');
          
          nQ = vv.N.(['Qg' stringIdx]);
          iQgVar = vv.i1.(['Qg' stringIdx]):vv.iN.(['Qg' stringIdx]);
          
          iBeta = vv.i1.(['Beta' stringIdx]):vv.iN.(['Beta' stringIdx]);
          if optns.gen.optimizeBaseP
              iPBase = vv.i1.('Pg'):vv.iN.('Pg');
          end
          
          %% get Pwind 
          wind = cs.wind(1+(i-1)*nCurtail:i*nCurtail)/baseMVA;
          
          % loop over all generators
          gIdx = 1;
          PVarCount = 1;
          PCurCount = 1;
          for j=1:ng
              if idxActive(j) && idxPQFactor(j) % generator active, add constraint
                  
                  qidx = iQgVar(gIdx);
                  % find corresponding active generation
                  if idxPVar(j) % variable P
                      pidx = iPgVar(PVarCount);
                      % add coefficients to A
                      A(Aidx,qidx) = 1;
                      A(Aidx,pidx) = -gen2(j,PQ_FACTOR);
                      vl(Aidx) = -inf; vu(Aidx) = 0;
                      Aidx = Aidx + 1;
                      A(Aidx,qidx) = -1;
                      A(Aidx,pidx) = -gen2(j,PQ_FACTOR);
                      vl(Aidx) = -inf; vu(Aidx) = 0;
                      Aidx = Aidx + 1;
                      
                      PVarCount = PVarCount + 1;
                  elseif idxPFix(j) % fixed P
                      % only constraint if P base are optimization
                      % variables
                      if optns.gen.optimizeBaseP
                          pidx = iPBase(j);
                          
                          A(Aidx,qidx) = 1;
                          A(Aidx,pidx) = -gen2(j,PQ_FACTOR);
                          vl(Aidx) = -inf; vu(Aidx) = 0;
                          Aidx = Aidx + 1;
                          A(Aidx,qidx) = -1;
                          A(Aidx,pidx) = -gen2(j,PQ_FACTOR);
                          vl(Aidx) = -inf; vu(Aidx) = 0;
                          Aidx = Aidx + 1;
                      end
                      
                  else % curtailable P
                      pidx = iBeta(PCurCount);
                      
                      A(Aidx,qidx) = 1;
                      A(Aidx,pidx) = gen2(j,PQ_FACTOR)*wind(PCurCount);
                      vl(Aidx) = -inf; vu(Aidx) = gen2(j,PQ_FACTOR)*wind(PCurCount);
                      Aidx = Aidx + 1;
                      A(Aidx,qidx) = -1;
                      A(Aidx,pidx) = gen2(j,PQ_FACTOR)*wind(PCurCount);
                      vl(Aidx) = -inf; vu(Aidx) = gen2(j,PQ_FACTOR)*wind(PCurCount);
                      Aidx = Aidx + 1;
                      
                      PCurCount = PCurCount + 1;     
                  end
                  
                  gIdx = gIdx + 1;
              end
          end
          
          
          %           % indices for PQcapability constraints
          %           if ~optns.gen.optimizeBaseP
          %               % add some constraints as variable limits
          %               idxPQcap = and( and( or(idxPVar, idxCurtail), idxActive ) , idxPQFactor);
          %               nPQcap = sum(idxPQcap);
          %           else
          %               % add all constraints as variable limits
          %               idxPQcap = and( idxActive, idxPQFactor );
          %               nPQcap = sum(idxPQcap);
          %           end
          
          %   OM = add_constraints(OM, NAME, A, L, U, VARSETS);
          vl = vl.';
          vu = vu.';
          
          om = add_constraints(om,['PQcap' stringIdx],A,vl,vu);
      end
  end
