function optns = check_opf_options(mpc,optns)



if ~isempty(optns.outputFile)   
    optns.fileID = fopen(optns.outputFile,'w');
else
    optns.fileID = '';
end

if size(optns.gen.windScenarios,1) > 1 && not(optns.gen.optimizeBaseP)
    if length(optns.gen.curtailableP) ~= size(optns.gen.windScenarios,1)
        error('The size of windScenarios does not match number of curtailable generators');
    end
end

if ~isempty(optns.gen.windProbabilities)
    assert(length(optns.gen.windScenarios)==length(optns.gen.windProbabilities),'Must specify probabilities for all wind scenarios');
end
%optns.mpopt.opf.violation = optns.foptions.TolCon;


% find indices of lines in transfer corridors
optns.corridorIdxs = cell(size(optns.transferCorridors));
for i=1:length(optns.transferCorridors)
    thisCorridor = optns.transferCorridors{i};
    optns.corridorIdxs{i} = line_idx(thisCorridor(:,1), ...
                               thisCorridor(:,2), mpc);
end

