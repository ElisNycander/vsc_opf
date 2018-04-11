function p = parse_penetration( casename )
% Parse the name of the case to determine the penetration level

[s,e] = regexp(casename, 'P0\d*\W','start','end');

if length(s) ~= 1
    warning('Could not find penetration level within ''casename''');
    p = [];
else
    pstring = casename(s+1:e-1);
end

p = str2double( pstring )/10^(length(pstring)-1);

    
