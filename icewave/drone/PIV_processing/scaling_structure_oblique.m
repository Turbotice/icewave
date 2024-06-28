function m = scaling_structure_oblique(S,x_pix,y_pix,t,x0,y0,H,alpha_0,focale,fps,Dt)
    
    % create a copy of the input structure 
    m = S;
    % Interpolate Vy for pixels that are at the center of each boxes, for
    % each frame
    Fy = griddedInterpolant({x_pix,y_pix,t},s.Vy);
    % compute Interpolant for Vz 
    [Fz] = backward_projection(Fy,y0,H,alpha_0,focale,fps,Dt); 

    [X_pix,Y_pix,T] = meshgrid(x_pix,y_pix,t);
    % Compute Vz 
    Vz = Fz(X_pix,Y_pix,T);
    Vz = permute(Vz,[2,1,3]);
    Vz = - Vz;
    m.Vz = Vz;
    
    % Compute values of x and y coordinates in real frame work
    [xreal,yreal] = projection_real_space(x_pix,y_pix,x0,y0,H,alpha_0,focale);
    m.x = xreal;
    m.y = yreal;
    m.t = (S.t - S.t(1)) * ft ; % convert to second 

    m.units.Vx = 'pix/frame';
    m.units.Vy = 'pix/frame';
    m.units.Vz = 'm/s';
    m.units.x = 'm';
    m.units.y = 'm';
    m.units.t = 's';
    
    m.xref = 'sensor center x0';
    m.yref = 'sensor center y0';
    
end 