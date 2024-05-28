%% Scripts enables to understand how the attenuation of a sinusoidal signal 
% modifies the Fourier spectrum 


%% Creation of a 1D decaying sinus
facq_x = 10;
x = (0:1/facq_x:100); 
y = (0:1/facq_x:80);
alpha = 0.1;
lambda = 5; 
phi = pi/4; % initial phase in rad
s_1D = exp(-alpha*x).*sin(2*pi*x/lambda + phi); 

figure(1)
plot(x,s_1D)

%% Create a 2D signal, sinusoidal in one direction 

s_2D = repmat(s_1D,size(y'))';
[Y,X] = meshgrid(y,x);

figure(2)
surf(X,Y,s_2D)
shading interp
view(2)
colormap(redblue)

%% Perfom 2D fft of this signal 
padding_bool = 1;
add_pow2 = 1;
[shifted_fft,fft_2D,kx,ky] = spatial_FFT(s_2D,padding_bool,add_pow2,facq_x);

[KY,KX] = meshgrid(ky,kx);

figure(3)
pcolor(kx,ky,abs(shifted_fft)')
shading interp
% axis([-2.5 2.5, -1 1])
% set(gca,'YScale','log')
% set(gca,'XScale','log')
% set(gca,'ColorScale','log')
colorbar()

%% Detect peaks 

zmax = max(abs(shifted_fft),[],'all');
norm_fft = abs(shifted_fft)/zmax;

% binarize the FFT spectrum
threshold = 0.8;
bin_fft = norm_fft;
bin_fft(norm_fft > threshold) = 1;
bin_fft(norm_fft <= threshold) = 0;

CC = bwconncomp(bin_fft');
stats = regionprops("table",CC,'Centroid');

figure(5) 
imagesc(bin_fft')
hold on 
plot(stats.Centroid(:,1),stats.Centroid(:,2),'or')

%% Fit a Lorentzienne curve on a peak 

i = 1;

row = round(stats.Centroid(i,1));
col = round(stats.Centroid(i,2));
% select a profile along kx 

window = 30;
profile_kx = abs(shifted_fft(row - window : row + window,col)).^2;
kx_tofit = kx(row - window : row + window);

[y_lorentz,p,~] = lorentzfit(kx_tofit,profile_kx',[],[],'3');

figure(6)
plot(kx_tofit,profile_kx,'.')
hold on 
plot(kx_tofit,y_lorentz,'r-')

alpha = p(3);
disp(['Alpha = ' num2str(alpha)])