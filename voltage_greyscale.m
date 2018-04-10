function voltage_greyscale(PH1,optns,table,color_factor,colors)

if optns.plotGreyscale
    
    for i=1:length(PH1)
        % Edge: Black - variable, Red - curtailable, None - Fixed
        % determine type of generator
        
        % determine location of generator
        % North: Blue, Central: Green, South: Yellow, External: Grey
        bus = table2array( table.Vm( i,'BUS' ) );
        bus2d = rem( bus,100 );
        if bus2d < 40 % North
            PH1(i).FaceColor = [1 1 1]*color_factor(1);
        elseif bus2d < 60 % Central
            PH1(i).FaceColor = [1 1 1]*color_factor(2);
        elseif bus2d < 70 % South
            PH1(i).FaceColor = [1 1 1]*color_factor(3);
        else % External
            PH1(i).FaceColor = [1 1 1]*color_factor(4);
        end
        %PH1(i).EdgeColor = [1 1 1];
    end
    
else
    
    for i=1:length(PH1)
        % Edge: Black - variable, Red - curtailable, None - Fixed
        % determine type of generator
        
        % determine location of generator
        % North: Blue, Central: Green, South: Yellow, External: Grey
        bus = table2array( table.Vm( i,'BUS' ) );
        bus2d = rem( bus,100 );
        if bus2d < 40 % North
            PH1(i).FaceColor = colors(1,:)*color_factor(1);
        elseif bus2d < 60 % Central
            PH1(i).FaceColor = colors(2,:)*color_factor(2);
        elseif bus2d < 70 % South
            PH1(i).FaceColor = colors(3,:)*color_factor(3);
        else % External
            PH1(i).FaceColor = colors(4,:)*color_factor(4);
        end
        %PH1(i).EdgeColor = [1 1 1];
    end

    
end
