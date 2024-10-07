function [Fp] = backward_projection(fun,y_0,h,alpha_0,f,fps,Dt)

    % Computes vertical velocity field using backward projection, 
    % from a pixel velocity field Vy(x,y,t) observed on a camera. 
    % This function takes as inputs : 
    % - fun : a function of coordinates : (x,y,t) in pixels and frames
    % - y_0 : pixel index of the central point of the camera sensor
    % - h : drone altitude in meter (above sea level)
    % - alpha_0 : inclination angle of the camera, angle to the horizontal 
    % - f : camera focal length
    % - fps : camera frame rate
    % - Dt : step between two images that are compared during PIV
    
    Fp = @(x,y,t)h*f*fun(x,y,t)./((y - y_0)*cos(alpha_0) + f*sin(alpha_0))...
    ./((y - y_0 + fun(x,y,t))*sin(alpha_0) - f*cos(alpha_0))*fps/Dt;
    
end 