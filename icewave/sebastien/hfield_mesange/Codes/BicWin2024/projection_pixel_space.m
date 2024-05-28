function [x,y] = projection_pixel_space(xp,yp,x_0,y_0,h,alpha_0,f)

    % Definition of X and Y in real framework
    xp = xp;
    yp = -yp;
    
    y = yp*f*sin(alpha_0)./(h/sin(alpha_0) - yp.*cos(alpha_0)) + y_0;
    x= xp/h.*(f*sin(alpha_0) + (y - y_0)*cos(alpha_0))+x_0;
    

end 