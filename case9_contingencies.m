function c = case9_contingencies


%% contingencies
% type: 0 - no contingency, 1- line trip,  2 - generator trip
% CONT_TYPE   CONT_IDX  LOAD_MULT LINE_MULT  CONT_PROB  
c = [
%    0    0      1       0                       % neutral scenario 
%	0     0    1.1      0                       % 10% load increase
%     0     0    0.9     0 % 10 % load decrease
    1     2    1       0   % trip line 2
% 	1     3    1       0   % trip line 3
% 	1     5    1       0   % trip line 5
% 	1     6    1       0   % trip line 6
% 	1     8    1       0   % trip line 8
% 	1     9    1       0   % trip line 9
	
% %     1     2    1       0.5   % weaken line 2    
     2     1    1       0   % trip generator 1
% 	2     2    1       0   % trip generator 1
%	2     3    1       0   % trip generator 1

%     2   1   1   0   1 % wind scenario nr 1, gen 1
%     2   1   1   0   2 % wind scenario nr 2, gen 1
%     1   2   1   0   1 % wind scenario nr 1, line 2
%     1   2   1   0   2 % wind scenario nr 2, line 2    
 %   2   1   1   0   0 % trip generator 1
];
%c = [];
