function [Fp] = backward_projection(fun,y_0,h,alpha_0,f,fps,Dt)



    %y_pix = repmat([1:1:length(yboite)]'*W/2 + 1,1,size(Vy,1))';
%     vx = Vx(:,:,i0);
    
%     y_prime = y_pix + vy;
%     x_prime = x_pix + vx;
    
    Fp = @(x,y,t)h*f*fun(x,y,t)./((y - y_0)*cos(alpha_0) + f*sin(alpha_0))...
    ./((y - y_0 + fun(x,y,t))*sin(alpha_0) - f*cos(alpha_0))*fps/Dt;
    
end 