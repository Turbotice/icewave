function [x,y] = projection_pixel_space(xreal,yreal,x_0,y_0,h,alpha_0,f)

    % Definition of x and y in pixel framework

    % Inputs : 
    % - xreal: array of x-coordinates in metric system
    % - yreal: array of y-coordinates in metric system
    % - x_0 : x-coordinate of camera sensor center (pixel system)
    % - y_0 : y-coordinate of camera sensor center (pixel system)
    % - h : drone altitude in meter (above sea level)
    % - alpha_0 : inclination angle of the camera, angle to the horizontal 
    % - f : camera focal length

    xreal = xreal;
    yreal = -yreal;
    
    y = yreal*f*sin(alpha_0)./(h/sin(alpha_0) - yreal.*cos(alpha_0)) + y_0;
    x = xreal/h.*(f*sin(alpha_0) + (y - y_0)*cos(alpha_0))+x_0;
    

end 