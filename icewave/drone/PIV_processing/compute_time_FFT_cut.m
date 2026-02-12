function [TF_cut_t] = compute_time_FFT_cut(V,nb_cuts,addpad)

% This function computes the smooth - time Fourier Transform of a
% velocity field. To do so, this function takes as arguments : 
% - V : the velocity matrix which is 3D [nx,ny,nt]
% - nb_cuts : the number of binin which we want to cut the whole video
% - addpad : the additional power of 2 used to compute the FFT on each bin
% the whole length of each sub_video will be : 2^(nextpow2(nt)  + addpad)

% The function returns the normalized, (non-scaled !!) matrix TF_cut_t. 
% This matrix has dimensions : [nx,ny,padding_length,nb_cuts] where the :
% - 3rd dimension is the frequency space computed over a given sub_video
% - 4th dimension is the sub_video index

%% FFT temporal smooth 
% Here we compute the temporal Fourier transform using cuts in time
% We cut the signal in nb_cuts bin along time, do the FFT_t for each bin
% and then compute an average of each FFT over all bins

% create cuts in the matrix, along time axis
[nx,ny,nt] = size(V);
r = rem(nt,nb_cuts); % reminder
q = (nt - r)/nb_cuts; % quotient 

% Cut the video matrix into nb_cuts sub_video matrices
% Each matrix has dimension over time of q
    for j = 1 : nb_cuts 
        V_cut(:,:,:,j) = V(:,:,(j-1)*q + 1 : j*q);
    end
    
% else % last cut is sufficiently big to be considered
%     nb_cuts = nb_cuts + 1;
%     for j = 1 : nb_cuts
%         Vx_cut(:,:,:,j) = Vx_s(:,:,(j-1)*q + 1 : j*q);
%     end 

% perform FFT temporal of each sub_video matrix 
    for j = 1:nb_cuts
        disp(j);
        signal_length = size(V_cut(:,:,:,j),3);
        padding_length = 2^(nextpow2(signal_length)  + addpad); % padding used to compute FFT_t
        disp(padding_length);

        TF_cut_t(:,:,:,j) = fft(V_cut(:,:,:,j),padding_length,3)/signal_length;
        % TF_cut_t(:,:,:,j) is the temporal Fourier Transform of the sub_video
        % number #j
    end 

end