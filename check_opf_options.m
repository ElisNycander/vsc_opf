function optns = check_opf_options(optns)


assert(isempty(setdiff(optns.gen.maxPg,optns.gen.fixPg)), ...
    'Error: ''maxPg'' must be subset of ''fixPg''' );

if ~isempty(optns.outputFile)   
    optns.fileID = fopen(optns.outputFile,'w');
else
    optns.fileID = '';
end
%optns.mpopt.opf.violation = optns.foptions.TolCon;


