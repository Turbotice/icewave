function [xreal,yreal] = projection_real_space(x,y,x_0,y_0,h,alpha_0,f)

    % Definition of x and y in real framework, camera sensor center is
    % taken as a reference 
    % Inputs : 
    % - x: array of x-coordinates in pixels
    % - y: array of y-coordinates in pixels
    % - x_0 : x-coordinate of camera sensor center
    % - y_0 : y-coordinate of camera sensor center
    % - h : drone altitude in meter (above sea level)
    % - alpha_0 : inclination angle of the camera, angle to the horizontal 
    % - f : camera focal length
    
    yreal = (y - y_0)*h/sin(alpha_0)./(f*sin(alpha_0) + (y - y_0)*cos(alpha_0));
    xreal = (x - x_0)*h./(f*sin(alpha_0) + (y - y_0)*cos(alpha_0));

    yreal = -yreal;

end 