
function m = scaling_structure(S,scale_V,fx,ft)
% This function scales data from a given structure, it also create a UTC
% time array if the structure presents an initial UTC time : 't0_UTC'

    m = S; % create a copy of the structure
    m.Vx = S.Vx * scale_V; % Convert to m/s
    m.Vy = S.Vy * scale_V; % Convert to m/s

    m.x = (S.x - S.x(1)) * fx ; % convert to meter
    m.y = (S.y - S.y(1)) * fx ; % convert to meter
    m.t = (S.t - S.t(1)) * ft ; % convert to second 

    m.units.Vx = 'm/s';
    m.units.Vy = 'm/s';
    m.units.x = 'm';
    m.units.y = 'm';
    m.units.t = 's';
    
    m.xref = 'bottom_left';
    m.yref = 'top_left';
    
    if isfield(S,'t0_UTC')
        time_s = milliseconds(m.t*1000); 
        t0 = S.t0_UTC;
        t0.Format = 'HH:mm:ss.SSS'; % displayed format
        m.UTC_t = t0 + time_s; % array of UTC time 
    else
        disp('Initial UTC time not found')
    end 

end 