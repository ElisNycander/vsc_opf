
function [f, df, d2f] = vscopf_f_minCurtail(x, om)
% Objective function: maximize sum of active generation from generators
% with gen2(:,PMAXIMIZE) == 1
define_constants;

%% unpack data
mpc = get_mpc(om);
vv = get_idx(om);

%% problem dimensions
nxyz = length(x);           %% total number of control vars of all types

%% grab Pg & Qg
%Pg = x(vv.i1.Pg:vv.iN.Pg);  %% active generation in p.u.

%widx = find(mpc.gen2(:,PMAXIMIZE));

% maximize voltage at all buses for base case
idx = zeros(1,nxyz);
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
idx(vv.i1.Qg) = 1;
        
idx = logical(idx);
    

% F: negative sum of wind
f = -sum(x(idx));
%x(idx)
if nargout > 1
    
  %% gradient

    df = zeros(nxyz,1);
    df(idx) = -1;
    
  %%  Hessian 
  if nargout > 2
       d2f = sparse(zeros(nxyz));
  end
end
