function m = scaling_structure_oblique(S,x_pix,y_pix,t,x0,y0,H,alpha_0,focale,fps,Dt)
    
    % create a copy of the input structure 
    m = S;
    % Interpolate Vy for pixels that are at the center of each boxes, for
    % each frame
    Fy = griddedInterpolant({x_pix,y_pix,t},S.Vy);
    % compute Interpolant for Vz 
    [Fz] = backward_projection(Fy,y0,H,alpha_0,focale,fps,Dt); 

    [X_pix,Y_pix,T] = meshgrid(x_pix,y_pix,t);
    % Compute Vz 
    Vz = Fz(X_pix,Y_pix,T);
    Vz = permute(Vz,[2,1,3]);
    Vz = - Vz;
    m.Vz = Vz;
    
    % Compute values of x and y coordinates in real frame work
    [X_pix,Y_pix] = meshgrid(x_pix,y_pix);
    X_pix = permute(X_pix, [2,1]);
    Y_pix = permute(Y_pix, [2,1]);
    [Xreal,Yreal] = projection_real_space(X_pix,Y_pix,x0,y0,H,alpha_0,focale);
    m.X = Xreal;
    m.Y = Yreal;
    m.t = (S.t - S.t(1)) / fps ; % convert to second 

    m.units.Vx = 'pix/frame';
    m.units.Vy = 'pix/frame';
    m.units.Vz = 'm/s';
    m.units.X = 'm';
    m.units.Y = 'm';
    m.units.t = 's';
    
    m.xref = 'sensor center x0';
    m.yref = 'sensor center y0';
    
    if isfield(S,'t0_UTC')
        time_s = milliseconds(m.t*1000); 
        t0 = S.t0_UTC;
        t0.Format = 'HH:mm:ss.SSS'; % displayed format
        m.UTC_t = t0 + time_s; % array of UTC time 
    else
        disp('Initial UTC time not found')
    end 
    
end 