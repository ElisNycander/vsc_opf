close all;
clear;
define_constants;
vscopf_define_constants;


%% OPTIONS
casename = 'case4';


%%
load([casename '.mat']);


% don't save plots
optns.saveFigures = 0;

% vv = get_idx(om);
% tab.Qg
% 180/pi*x(vv.i1.Va1:vv.iN.Va1)
%restab.Vm
restab.Pg
%restab.Beta
%restab.Curtail
%restab.ExpCurtail
restab.Qg
%restab.S
%restab.Slam
%restab.lamInfo
restab.PQ


plot_results(restab,optns);