function [alpha,C,d] = exponential_coeff(A,x)
% Perform exponential fit of a profile
% Inputs : - A : 1D array, Amplitude
%          - x : 1D array, x coordinate
% 
% Outputs : - alpha : scalar, exponential coeff
%           - C : scalar, pre-factor
%           - d : scalar, distance to exponential fit 

    log_A = log10(A); % take the log10 of the amplitude of freq i
    A_fit = log_A; 
    p = polyfit(x,A_fit,1); % fit log_A by a degree 1 polynome
    alpha = log(10)*p(1); % get attenuation coefficient in m^-1
    C= 10.^p(2); % prefactor
    
    disp(['alpha = ' num2str(alpha) ' m-1'])
    y_poly = 10.^polyval(p,x);

    d = sum((y_poly - A).^2)/sum(A.^2); % distance to the fit
    disp(['Distance to fit : ' num2str(d) ''])


end 