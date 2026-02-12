function [E,f,k,shift] = space_time_spectrum_Efk(V,add_pow2,facq_t,facq_x,padding_time)

% Computes space-time amplitude spectrum of the 3D-array V
% Inputs: - V : 3D-array, dimensions (nx,ny,nt)
%         - add_pow2 : additional padding if needed, 1x3 array 
%         - facq_x : acquisition frequency in box/meter
%         - facq_t : acquisition frequency in frame/meter
%         - padding_time : boolean to pad time or not 
% Outputs : - E : space-time spectrum, order of dimensions (f,k)
%           - f : frequency 1D-array 
%           - k : wave vector 1D-array 
%           - shift : array of shifted 2D-FFT for each positive frequency
%           (kx,ky,f)

% FFT 3D of the velocity field 
N = size(V);
padding = 2.^(nextpow2(N) + add_pow2); % padding for each dimension
nopadding_time = ~padding_time;

if nopadding_time % if we do not want to pad time 
    padding(3) = N(3);
end 

disp('Computing FFT 3D')
FFT = fftn(V,padding)/numel(V); % FFT 
disp('FFT 3D computed')

kx = 2*pi*facq_x*(-padding(1)/2:padding(1)/2-1)/padding(1);
ky = -2*pi*facq_x*(-padding(2)/2:padding(2)/2-1)/padding(2);

% keep only positive frequencies 
FFT_positive = FFT(:,:,1:padding(3)/2 +1);
FFT_positive(:,:,2:end-1) = 2*FFT_positive(:,:,2:end-1);

f = facq_t*(0:padding(3)/2)/padding(3);

% FFT shift for all dimensions
shift = fftshift(fftshift(FFT_positive,2),1);
disp('FFT shifted')

% Radial average for each frequency

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

k = 2*pi*R_tics*facq_x/padding(1);

end 
