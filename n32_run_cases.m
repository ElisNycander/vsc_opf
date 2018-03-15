close all;
clear;

rerun_simulations = 1;
base_identifier = 'run1';

cases = {'B1','B2','C1','C2','D1','D2'};
PQ_constraints = [0 1];
Q_wind = [0 1];


%% run simulations
ijk = 1;
for i=1:length(cases)
    for j=1:length(PQ_constraints)
        for k=1:length(Q_wind)
            
            caseNames{ijk} = cases{i};
            if PQ_constraints(j) == 1
                caseNames{ijk} = [caseNames{ijk} '_PQ'];
            end
            if Q_win(k) == 1
                caseNames{ijk} = [caseNames{ijk} '_Qwind'];
            end
            % add 
            caseNames{ijk} = [base_identifier '_' caseNames{ijk}];
            
            if rerun_simulations
                run_opf_n32_fcn(caseNames{ijk},cases{i},PQ_constraints(j),Q_wind(k));
            end
            
            ijk = ijk + 1;
        end
    end 
end

%% collect results
for i=1:length(caseNames)
    
    
    
end
