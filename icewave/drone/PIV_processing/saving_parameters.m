function saving_parameters(m,filename,parameters_list)
% Save relevant parameters associated to a PIV movie in a file 

    fileID = fopen(filename,'a');
    
    for k = 1 : length(parameters_list)
        param = parameters_list{k};
        if isfield(m.units,param)
            unit = getfield(m.units,param);
            %     unit = units_list{k};
            formatSpec = [param '  %4.3f  ' unit ' \n'];
            fprintf(fileID,formatSpec,m.(param));
        else 
            disp( 'No unit found for this variable')
            formatSpec = [param '  %4.3f  \n'];
            fprintf(fileID,formatSpec,m.(param));
        end 
            
    end 

    fclose('all');
end 

