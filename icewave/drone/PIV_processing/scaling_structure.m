
function m = scaling_structure(S,scale_V,fx,ft)
% This function scales data from a given structure, it also create a UTC
% time array if the structure presents an initial UTC time : 't0_UTC'

    m = S; % create a copy of the structure
    m.Vx = S.Vx * scale_V; % Convert to m/s
    m.Vy = S.Vy * scale_V; % Convert to m/s

    m.x = S.x * fx ; % convert to meter
    m.y = S.y * fx ; % convert to meter
    m.t = (S.t - S.t(1)) * ft ; % convert to second 
    
    [X,Y] = meshgrid(m.x,m.y);
    m.X = permute(X,[2,1]);
    m.Y = permute(Y,[2,1]);
    
    m.units.Vx = 'm/s';
    m.units.Vy = 'm/s';
    m.units.x = 'm';
    m.units.y = 'm';
    m.units.X = 'm';
    m.units.Y = 'm';
    m.units.t = 's';
    
    m.xref = 'bottom_left_box_center';
    m.yref = 'top_left_box_center';
    
    if isfield(S,'t0_UTC')
        time_s = milliseconds(m.t*1000); 
        t0 = S.t0_UTC;
        t0.Format = 'HH:mm:ss.SSS'; % displayed format
        m.UTC_t = t0 + time_s; % array of UTC time 
    else
        disp('Initial UTC time not found')
    end 
    
    % convert UTC_t array in array of character 
    m.UTC_t = char(m.UTC_t);


end 