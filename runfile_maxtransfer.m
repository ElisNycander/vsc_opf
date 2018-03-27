clear;
close all;
define_constants;

mpc = Nordic32;

buses = [
    4031
    4032
    4022
    4021
    4012
    4011
    4041
    4042
    4044
    4043
    4045
    4046
    ];

nbus = length(buses);

buses2d = rem(buses,100);
idxBusNorth = buses2d < 40;
idxBusSouth = buses2d > 40;

Qvalues = 0:100:1e3;

%cs = struct();

max_transfer = [];
curve_labels = cell(size(buses));
base_transfer = [];
base_Q = [];

for i=1:nbus
    
    cs(i) = maximum_transfer_vs_qlim(mpc,buses(i),Qvalues,'Nordic 32 original');
    
    max_transfer = [max_transfer; cs(i).maxTransfer];
    base_transfer = [base_transfer; cs(i).baseTransfer];
    base_Q = [base_Q cs(i).baseQ];
    curve_labels{i} = ['Bus ' num2str(buses(i))];
end
%curve_labels{nbus+2} = 'base cases';
%curve_labels{nbus+2} = 'current limit';


figure; 
PH1 = plot(Qvalues,max_transfer); hold on; grid on;
plot(base_Q,base_transfer,'ko');

xlabel('Q_{max} (MW)');
ylabel('Transfer (MW)');
LH1 = line(xlim, [cs(1).baseTransfer cs(1).baseTransfer]);
LH1.Color = [0 0 0];
LH1.LineStyle = '--';
title('Maximum transfer capacity vs Q_{max} at different buses');
legend(PH1, curve_labels);

for i=1:length(PH1)
    
    if idxBusNorth(i) == 1
        PH1(i).LineWidth = 2;
    end
    
end

% Plot sensitivity of transfer to increase in Qmax

dTxdQ = zeros(nbus,1); % MW/MW

for i=1:nbus
    % find point after base transfer
    idxQ1 = find( Qvalues > base_Q(i) , 1 );
    dTxdQ(i) = ( max_transfer( i,idxQ1 )- base_transfer( i ) ) / ...
        ( Qvalues( idxQ1 ) - base_Q( i ) );
end

figure;
bar( categorical(curve_labels), dTxdQ.' );
ylabel('d(Transfer capacity)/dQ (MW/MW)');
title('Sensitivity of transfer capacity to additional Q capability');


