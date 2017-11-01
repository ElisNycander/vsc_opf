function c = case9_contingencies


%% contingencies
% type: 0 - no contingency, 1- line trip,  2 - generator trip
% CONT_TYPE   CONT_IDX  LOAD_MULT LINE_MULT   
c = [0     0    1.1     0 % 10% load increase
    1     2    1       0   % trip line 2
     1     2    1       0.5   % weaken line 2    
    2     1    1       0   % trip generator 1
];
%c = [];
