function [E,f,k,shift] = get_A_fk(V,add_pow2,facq_t,facq_x)
    % This function computes the amplitude A in the Fourier space (f,k) of a
    % matrix V, size [Nx,Ny,Nt]. It performs an FFT along each direction (2
    % directions of space + 1 direction of time), than for each frequency,
    % performs an average over orientations of wave vector k. 

    % It takes as arguments : 
    % - field V, matrix of dimension [Nx,Ny,Nt]
    % - add_pow2, additional padding for each dimension, array  1 x 3
    % - facq_t, acquisition frequency for time dimension 
    % - facq_x, acquisition frequency for space dimensions

    % FFT 3D of the velocity field 
    N = size(V);
    padding = 2.^(nextpow2(N) + add_pow2); % padding for each dimension

    padding(1) = max(padding(1),padding(2)); % keep only biggest padding over space
    padding(2) = max(padding(1),padding(2));

    disp('Computing FFT over each dimension')
    FFT = fftn(V,padding)/numel(V); % FFT scaled using initial number of elements in V
    disp('FFT computed')

    %%
    kx = 2*pi*facq_x*(-padding(1)/2:padding(1)/2-1)/padding(1);
    ky = -2*pi*facq_x*(-padding(2)/2:padding(2)/2-1)/padding(2);

    %% keep only positive frequencies 
    FFT_positive = FFT(:,:,1:padding(3)/2 +1);
    FFT_positive(:,:,2:end-1) = 2*FFT_positive(:,:,2:end-1);

    f = facq_t*(0:padding(3)/2)/padding(3);

    %% FFT shift for dimensions in space 
    shift = fftshift(fftshift(FFT_positive,2),1);
    disp('FFT shifted')

    %% Radial average for each frequency

    radial_step = 1;
    % center of the 2D-FFT
    x0 = padding(1)/2 + 1;
    y0 = padding(2)/2 + 1;

    %loop over all frequencies 
    for i0 = 1:size(shift,3)
        % demodulated 2D-FFT
        img = shift(:,:,i0);

        [R_tics, R_profile] = radialavg2(abs(img),radial_step,x0,y0);
        E(i0,:) = R_profile; % Amplitude for value (f,k)
    end 

    disp('Amplitude computed')
    % wave vector array, size : padding(1)*sqrt(2)/2
    k = 2*pi*R_tics*facq_x/padding(1); 

    disp('DONE.')
end 