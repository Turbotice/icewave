
function [converted_s] = conversion_structure_units(s,keys,initial_units,IS_units,coef_conversion)
% Convert field units of a given structure 
% Inputs : - s : structure containing fields and a field 'units' with the
% same fields that appear in s
%          - keys : cell arrray of strings, fields to be converted
%          - initial_units : cell array of strings, initial unit of the
%          corresponding key 
%          - IS_units : cell array of strings, final unit of the
%          corresponding key
%          - coef_conversion : array of double, conversion coefficients
%          for each key, from initial_units to IS_units. 
%
% Outputs : - converted_s : structure containing converted fields 

    if ~isfield(s,'units')
        disp(' The studied structure does not contain any field "units" ')
    end

    for key_idx = 1 : length(keys)
        key = keys{key_idx};
        
        if isfield(s,key) % key is a field of the structure
    
            if strcmp(s.('units').(key),initial_units{key_idx})
                s.(key) = s.(key)*coef_conversion(key_idx);
                s.('units').(key) = IS_units{key_idx};
            end 
    
        else 
            disp([key ' is not a field of the studied structure'])
        end 
    
    end 

    converted_s = s;
end