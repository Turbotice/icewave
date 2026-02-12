function [omegaN] = bound_harmonicN(k,N,h_w)
    % Computes the pulsation omegaN for a given k, along the dispersion
    % relation of the harmonic N
    
    omegaN = sqrt(9.81*N.*k.*tanh(h_w.*k./N));
end