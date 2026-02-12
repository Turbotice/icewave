function [shifted_fft,fft_2D,ky,kx] = spatial_FFT(img,padding_bool,add_pow2,fx)
% This function computes the 2D-FFT of an image
% It takes as arguments: 
% - the image : img, size : (Ny,Nx)
% - padding_bool : a boolean variable to know if we use padding for the
% computation of the 2D-FFT or not
% - add_pow2 : The additional power of 2 to pad the array with more zeros
% - fx : the scale used in pixel/meter

    [Ny, Nx] = size(img);
    Nb_elements = Nx*Ny; % number of elements of matrix a
    if padding_bool 
            % generate a 2D - hamming function
            ham_window = window2(Ny,Nx,@hamming);
            
            padding_x = 2^(nextpow2(Nx) + add_pow2);
            padding_y = 2^(nextpow2(Ny) + add_pow2);
            fft_2D = fft2(img.*ham_window,padding_y,padding_x);
            disp('Padding used')
    else 
            % generate a 2D - hamming function
            ham_window = window2(Ny,Nx,@hamming);
            fft_2D = fft2(img.*ham_window);
            padding_x = Nx;
            padding_y = Ny;
            disp('No padding')
    end

    fft_2D = fft_2D/Nb_elements; % normalization by the initial number of elements
    % #################################
    % SCALING : create kx and ky arrays
    % #################################
    kx = -2*pi*fx*(-padding_x/2:padding_x/2-1)/padding_x;
    ky = -2*pi*fx*(-padding_y/2:padding_y/2-1)/padding_y;

    shifted_fft = fftshift(fft_2D);

end 