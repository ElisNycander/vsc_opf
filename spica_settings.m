
%% lines that constitute bottleneck
spica_bottleneck = [
    % from      to
    4031        4041
    4032        4044
    4032        4042
    4021        4042
    ];


%% list with contingencies 
spica_contingencies = [
    % bus1      bus2    type (1- trip line, 2- trip generator)
    
%     % base case
%     0           0       0  
% 
%     % tripping of bottleneck lines
%      [spica_bottleneck ones(size(spica_bottleneck,1),1)];
% %     
% %     % tripping of other lines
%        4011        4021    1 
% %     
% %     % tripping of generators
%     1042        0       2
%     1043        0       2              
%     4042        0       2
%    4047        0       2
  4051        0       2
    
];

%%
% buses where generation should be increased when solving cpf
% note: relative distribution of increase can be according to specified +PG
% weights or the same as the initial generation distribution among buses
spica_gen_increase = [
    % bus   +PG (MW)
    1012    11
    1013    11
    1014    11
    1021    11
    1022    11
    2032    11
    4011    11
    4012    11
    4021    11
    4031    11
    ];

% buses where load is increased for cpf
spica_load_increase = [
    % bus +PD (MW)
    1041    10
    1042    10
    1043    10
    1044    10
    1045    10
    41      10
    42      10
    43      10
    46      10
    47      10
    51      10
];

    
