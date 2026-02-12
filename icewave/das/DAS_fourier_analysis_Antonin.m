%% Data reading
clear all
f=rdir('*.h5');
%h5disp(f{1})
raw_data=h5read(f{4},"/Acquisition/Raw[0]/RawData");
%%
% clean up of the data.
% Extraction of the coupled part of the fiber. From spatial pixel 18+5 to 18+33
% Subsampling from 40 kHz to 1 kHz
data=raw_data(18+5:18+33,1:40:end);
imagesc(data)
%% Fourier analysis
%
Spatial_padding=128;
time_length=size(data,2)/10;
Y=fftshift(fft2(data(:,1:10:end)-mean(data,2),Spatial_padding,size(data,2)/10)); % Subsampling to 100 Hz
%imagesc(abs(Y(:,size(data,2)/2+1-100:size(data,2)/2+1+100)))
% Scales in physical units
% temporal sampling : 1kHz;
%f_s=1000;
f_s=100; % Subsampling above
Total_length=time_length;
total_duration=time_length/f_s;
delta_f=1/total_duration;
f_axis= ([1:time_length]-time_length/2-1)*delta_f;
% Spatial sampling : 4.785714 m/pix
Spatial_scale=4.785714;
Total_spatial_length = Spatial_padding *Spatial_scale % Includes padding
delta_k=1/Total_spatial_length;
k_axis=([1:Spatial_padding]-Spatial_padding/2-1)*delta_k*2*pi;
%% Plot de la fft2
imagesc(-k_axis,f_axis,abs(Y'))