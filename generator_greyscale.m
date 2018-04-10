function generator_greyscale(PH1,optns,table,color_factor,colors)


for i=1:length(PH1)
    % Edge: Black - variable, Red - curtailable, None - Fixed
    % determine type of generator
    if ismember( i,optns.gen.curtailableP ) % curtailable
        PH1(i).EdgeColor = [0 0 0];
        PH1(i).LineWidth = 0.8;
    elseif ismember( i,optns.gen.fixedP ) % fixed
        PH1(i).EdgeColor = [1 1 1];
    else % variable
        PH1(i).EdgeColor = [0 0 0];
        PH1(i).LineWidth = 0.8;
    end
    
    
    
    % determine location of generator
    % North: Blue, Central: Green, South: Yellow, External: Grey
    genbus = table2array( table.Qg( i,'BUS' ) );
    genbus2d = rem( genbus,100 );
    if optns.plotGreyscale
        if genbus2d < 40 % North
            PH1(i).FaceColor = [1 1 1]*color_factor(1);
        elseif genbus2d < 60 % Central
            PH1(i).FaceColor = [1 1 1]*color_factor(2);
        elseif genbus2d < 70 % South
            PH1(i).FaceColor = [1 1 1]*color_factor(3);
        else % External
            PH1(i).FaceColor = [1 1 1]*color_factor(4);
        end
        
    else
        if genbus2d < 40 % North
            PH1(i).FaceColor = colors(1,:)*color_factor(1);
        elseif genbus2d < 60 % Central
            PH1(i).FaceColor = colors(2,:)*color_factor(2);
        elseif genbus2d < 70 % South
            PH1(i).FaceColor = colors(3,:)*color_factor(3);
        else % External
            PH1(i).FaceColor = colors(4,:)*color_factor(4);
        end
        
    end
end
