function double_x_axis(FH,oldaxis,newaxis)
% Add extra labels to axis, assuming linear conversion from old scale to
% new scale. 
% Inputs: 
% FH - figure handle
% oldaxis - old scale 
% newscale - new scale

figure(FH);
%% for script
% figure(12);
% oldaxis = penetration;
% newaxis = windscale;

%% function
ymin = min(ylim);
yrng = diff(ylim);
xt = get(gca,'XTick');

% convert to column vectors
if size(oldaxis,1) > 1
    if size(oldaxis,2) > 1
        error('old axis not a vector');
    else
        oldaxis = oldaxis.';
    end
end
if size(newaxis,1) > 1
    if size(oldaxis,2) > 1
        error('new axis not a vector');
    else
        newaxis = newaxis.';
    end
end
if length(newaxis) ~= length(oldaxis)
    error('axes have different lengths')
end

% extend mapping outside given interval
oldaxis = [2*oldaxis(1)-oldaxis(end) oldaxis 2*oldaxis(end)-oldaxis(1)];
newaxis = [2*newaxis(1)-newaxis(end) newaxis 2*newaxis(end)-newaxis(1)];

% map ticks from old to new axis
xtn = interp1( oldaxis,newaxis,xt );

set(gca, 'XTickLabel',[])
nrxt = size(xt,2); 
latvct = linspace(xt(1), xt(end), nrxt);
lat_txt = cellstr(strsplit(num2str(latvct,'%.2f '),' '));
lat_txt{1} = ['%  ' lat_txt{1}];
text([(xt(1)-xt(2))*0.5 xt], ones(1,nrxt+1)*ymin-0.05*yrng, ['% ' lat_txt],'HorizontalAlignment','center')
lngvct = linspace(xtn(1), xtn(end), nrxt);
lng_txt = cellstr(strsplit(num2str(lngvct,'%.0f '),' '));
lng_txt{1} = ['MW ' lng_txt{1}];
text([(xt(1)-xt(2))*0.5 xt], ones(1,nrxt+1)*ymin-0.1*yrng, ['MW' lng_txt],'HorizontalAlignment','center')
