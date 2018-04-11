function change_markers(PH, marker_style, marker_size)

if length(PH) == length(marker_style) 
for i=1:length(PH)
    PH(i).Marker = marker_style{i};
    PH(i).MarkerSize = marker_size(i);
end

elseif rem(length(PH), length(marker_style) ) == 0
    
    nscenarios = length(PH) / length(marker_style);
    
    % take select colors from grayscale
    
    cfactor = linspace(0.7,0,nscenarios);
    
    ij = 1;
    for i=1:length(marker_style) % cases
        for j=1:nscenarios % scenarios
            PH(ij).Marker = marker_style{i};
            PH(ij).MarkerSize = marker_size(i);
            PH(ij).Color = [1 1 1]*cfactor(j);
            ij = ij + 1;
        end
    end
    
    
    
else
    warning('Number of plots is not multiple of number of cases given in marker_style');
end