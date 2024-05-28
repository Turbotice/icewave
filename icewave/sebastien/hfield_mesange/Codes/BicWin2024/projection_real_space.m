function [xp,yp] = projection_real_space(x,y,x_0,y_0,h,alpha_0,f)

    % Definition of X and Y in real framework
    
    yp = (y - y_0)*h/sin(alpha_0)./(f*sin(alpha_0) + (y - y_0)*cos(alpha_0));
    xp = (x - x_0)*h./(f*sin(alpha_0) + (y - y_0)*cos(alpha_0));
    
    xp = xp;
    yp = -yp;

end 