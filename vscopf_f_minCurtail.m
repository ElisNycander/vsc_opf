
function [f, df, d2f] = vscopf_f_minCurtail(x, om)
% Objective function: maximize sum of active generation from generators
% with gen2(:,PMAXIMIZE) == 1
define_constants;

%% unpack data
mpc = get_mpc(om);
vv = get_idx(om);
wind = mpc.contingencies.wind/mpc.baseMVA; % wind scenarios
prob = mpc.contingencies.probabilities;
nCurtail = length(find(mpc.gen2(mpc.gen2(:,PTYPE) == PCUR)));
%% problem dimensions
nxyz = length(x);           %% total number of control vars of all types

idx = false(1,nxyz);

for i=2:mpc.contingencies.N % note: no curtailment for base case
   sidx = num2str(i); 
   idx(vv.i1.(['Beta' sidx]):vv.iN.(['Beta' sidx])) = true; 
end
prob = zeros(nCurtail*(mpc.contingencies.N-1),1);
for i=1:mpc.contingencies.N-1
    prob(1+(i-1)*nCurtail:i*nCurtail) = mpc.contingencies.probabilities(i+1);
end
f = sum( x(idx) .* wind(nCurtail+1:end) .* prob );



% idx(vv.i1.Vm:vv.iN.Vm) = 1;
% idx = logical(idx);

% maximize injected active power
% if isfield(vv.i1,'Pg')
%     idx(vv.i1.Pg:vv.iN.Pg) = 1;
% end
% flag = 1;
% nr = 1;
% while flag
%     if isfield(vv.i1,['Pg' num2str(nr)])
%         idx(vv.i1.(['Pg' num2str(nr)]):idx.iN.(['Pg' num2str(nr)])) = 1;
%     else
%         flag = 0;
%     end
% end

% maximize Qg from generator 1
% idx(vv.i1.Qg) = 1;
%         
% idx = logical(idx);
%     
% 
% % F: negative sum of wind
% f = -sum(x(idx));
%x(idx)
if nargout > 1
    
  %% gradient

    df = zeros(nxyz,1);
    df(idx) = wind(nCurtail+1:end) .* prob;
    
  %%  Hessian 
  if nargout > 2
       d2f = sparse(zeros(nxyz));
  end
end
