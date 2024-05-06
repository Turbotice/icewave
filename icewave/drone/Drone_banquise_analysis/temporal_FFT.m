function [FFT_t,TF_spectrum,f] = temporal_FFT(H_ww,padding_bool,add_pow2,fps)
    % Description
    % This function computes the temporal FFT of the height matrix

    % It returns :
    % - FFT_t : a 3D matrix (nx,ny,nf)
    % - TF_spectrum : a 1D matrix, which represents the spacial average of
    % FFT_t 
    % - f : the frequency array [nf]

    % It takes as arguments :
    % - H_ww : the height matrix, obtained from FCD demodulation [nx,ny,nt]
    % - padding_bool : a boolean to know if we want to padd or not
    % - add_pow2 : additional power of 2 for padding
    % - fps : the frame rate used

    % FFT temporal spectrum of the image 
    original_length = size(H_ww,3); % original length of the signal
    padding_length = 2^(nextpow2(original_length) + add_pow2);
    
    % averaging over time 
    H_mean = mean(H_ww,3);
    H_ww = H_ww - H_mean;

    if padding_bool 
        FFT_t = fft(H_ww,padding_length,3);
        N = padding_length ;
        disp('Padding used')
    else 
        FFT_t = fft(H_ww,[],3);
        N = original_length;
        disp('No padding')
    end

    FFT_t = FFT_t/original_length; % normalization of the FFT
    TF_inter = squeeze(mean(mean(abs(FFT_t),2),1));

    TF_spectrum = TF_inter(1:N/2+1);
    TF_spectrum(2:end-1) = 2*TF_spectrum(2:end-1); % multiply by 2 for peaks that are both in positive an negative frequencies
    f = fps*(0:(N/2))/N; % frequency array
    
    % keep only one half of the spectrum
    FFT_t = FFT_t(:,:,1:N/2+1);
    FFT_t(:,:,2:end-1) = 2*FFT_t(:,:,2:end-1);

end 
