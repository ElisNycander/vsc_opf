function str = create_cpf_str(desc,mpcc,mpcb,mpc)
% Create string containing information about node conversion from cpf.
define_constants;

 str = sprintf('-------------------------------\n');
    if isfield(mpcc,'cpf')
        
        
    str = [str sprintf([desc '\n' mpcc.cpf.done_msg '\nBase transfer: %0.0f MW\nMax transfer: %0.0f MW\n'], ...
                        [mpcb.transfer(1) mpcc.transfer(1)])];
   
        str = [str sprintf('\nBase case node conversion:\n')];
        
        % find buses that are converted in base case
        ib = find(mpcb.bus(:,BUS_TYPE) ~= mpc.bus(:,BUS_TYPE));
        for ii=1:length(ib)
            bus = mpc.bus(ib(ii),1);
            if mpc.bus(ib(ii),BUS_TYPE) == PV && mpcb.bus(ib(ii),BUS_TYPE) == PQ
                str = [str sprintf('Bus %d from PV to PQ:\n',bus)];
            elseif mpc.bus(ib(ii),BUS_TYPE) == PV && mpcb.bus(ib(ii),BUS_TYPE) == 3
                str = [str sprintf('Bus %d from PV to SLACK:\n',bus)];
            elseif mpc.bus(ib(ii),BUS_TYPE) == 3 && mpcb.bus(ib(ii),BUS_TYPE) == PV
                str = [str sprintf('Bus %d from SLACK to PV:\n',bus)];
            elseif mpc.bus(ib(ii),BUS_TYPE) == 3 && mpcb.bus(ib(ii),BUS_TYPE) == PQ
                str = [str sprintf('Bus %d from SLACK to PQ:\n',bus)];
            else
                error('Unknown bus conversion at bus %d',bus);
            end
            % find generators at this bus
            gi = find(mpc.gen(:,GEN_BUS)==bus);
            pg = mpcb.gen(gi,QG);
            for iii=1:length(gi)
                str = [str sprintf('Generator %d: Q = %0.1f MW\n',[gi(iii) pg(iii)])];
            end
        end
        str = [str sprintf('\nCPF node conversion:\n')];
        
        events = mpcc.cpf.events;
        for ii=1:length(events)
            str = [str events(ii).msg sprintf('\n')];
        end
    else
        str = [str sprintf([mpcb.desc '\nUnsuccesful\n '])];
    end
    str = [str sprintf('-------------------------------')];