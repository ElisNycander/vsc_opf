function check_opf_options(optns)


assert(isempty(setdiff(optns.gen.maxPg,optns.gen.fixPg)), ...
    'Error: ''maxPg'' must be subset of ''fixPg''' );


%optns.mpopt.opf.violation = optns.foptions.TolCon;


