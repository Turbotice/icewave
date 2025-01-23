function [field_star,x_star,y_star,X_line,Y_line] = orientation_parallel2propagation(field,theta,fx,L0)
    % Select a rectangular window aligned with the propagation direction 
    % Within this window, create an interpolated field, and associated
    % coordinates : x_star (direction of propagation), y_star (oriented
    % upward) 
    % !!! WE CONSIDER ONLY WAVE FIELDS PROPAGATING ALONG + x-direction 
    % Inputs : 
    % - field : 2D array [Nx,Ny], can be complex values 
    % - theta : angle of the propagation direction, defined between
    % [-pi;pi], positive in counterclockwise 
    % - fx : scaling factor in meter/box 
    % - L0 : length of the window along the direction of propagation
    % 
    % Outputs : 
    % - field_star : 2D array interpolated [Nx,Ny] along the direction
    % given by angle theta
    % - x_star : curvilinear coordinate (x-coordinate in field_star
    % framework)
    % - y_star : y-coordinate in field_star framework 
    
    
    x = (1:1:size(field,1))*fx;
    y = (1:1:size(field,2))*fx;

    % Build several segments for theta < 0
    if theta < 0
        
        % initial line 
        x0 = x(1); % minimal value along x axis
        y0 = L0*sin(abs(theta)) + y(1);
        ds = fx; % step of curvilinear coordinate

        s = (0:ds:L0); % curvilinear coordinate

        Nb_lines = floor((y(end) - y0)/(ds*cos(theta))); % number of lines to draw 
        X_line = zeros(length(s),Nb_lines);
        Y_line = zeros(length(s),Nb_lines);

        for j = 1:Nb_lines
            x0 = x0 + ds*sin(abs(theta));
            y0 = y0 + ds*cos(theta);

            X_line(:,j) = x0 + s*cos(theta); % #1 : s-coordinate, #2 line index 
            Y_line(:,j) = y0 + s*sin(theta);
        end 

        % Interpolate 
        F = griddedInterpolant({x,y},field); % Interpolate demodulated field 
        % interpolate along new grid 
        field_star = F(X_line, Y_line);

        x_star = s; 
        y_star = (0:ds:(y(end) - (L0*sin(abs(theta)) + y(1)))/cos(theta)-ds);
    
    else % for theta > 0
        x0 = (y(end) - L0*sin(theta) - y(1))*tan(theta);
        y0 = y(1);
        ds = fx;
        
        s = (0:ds:L0); % curvilinear coordinate
        Nb_lines = floor((y(end) - L0*sin(theta) - y(1))/(ds*cos(theta))); % number of lines to draw 
        X_line = zeros(length(s),Nb_lines);
        Y_line = zeros(length(s),Nb_lines);
        
        for j = 1:Nb_lines
            x0 = x0 - ds*sin(abs(theta));
            y0 = y0 + ds*cos(theta);

            X_line(:,j) = x0 + s*cos(theta); % #1 : s-coordinate, #2 line index 
            Y_line(:,j) = y0 + s*sin(theta);
        end 

        % Interpolate 
        F = griddedInterpolant({x,y},field); % Interpolate demodulated field 
        % interpolate along new grid 
        field_star = F(X_line, Y_line);

        x_star = s; 
        y_star = (0:ds:(y(end) - L0*sin(abs(theta)) - y(1))/cos(theta) - ds);
        
    end
    
    
    
end 