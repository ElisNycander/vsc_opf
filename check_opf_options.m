function optns = check_opf_options(optns)



if ~isempty(optns.outputFile)   
    optns.fileID = fopen(optns.outputFile,'w');
else
    optns.fileID = '';
end
%optns.mpopt.opf.violation = optns.foptions.TolCon;


