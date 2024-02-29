%% Clear 
clear all;
close all;

%% Loading structure obtained after PIV processing and post-processing

base = '//192.168.1.70/Share/Data/0211/Drones/Fulmar/matData/';
%base = 'E:/PIVlab_drone/matdata/DJI_0308_Dt4_W64_full/';

filename = 'PIV_processed_i04500_Dt4_b1_W32_full_total_processed.mat';
matname = [base filename];
%%
disp('Loading Data..');
load(matname);
disp('Data loaded');
%% Creates a folder where to save the generated plots
base_fig = '//192.168.1.70/Share/Data/0211/Drones/Fulmar/matData/';
%base_fig = 'E:/PIVlab_drone/';
fig_folder = [base_fig 'SWO_FUL_20240211T203233UTC/Plots/'];
if ~exist(fig_folder)
    mkdir(fig_folder)
end

%% Scaling 
fe = 29.97; % Frame rate in Hz
nb_pass = m.s_param{6,2};
W = 36; % size of the window for DIC algorithm
font_size = 13;

% ##########################################
L_x = 3840; % size of the image in pixel, along x-axis
h_drone = 140; % height of the drone in meter
theta_x = 32.75; % semi AFOV of the drone, along x-axis, in °

fx_pix = L_x/(2*h_drone*tan(theta_x*pi/180)); % scale in pixels / meter
fx = fx_pix*2/W;
% ##########################################

% fx = 0.8857; % for W = 64; 
%fx = 0.8857*2; % spatial scale in boxes/meter
Dt = 4; % step between two frames that were compared during the PIV algorithm 

scale_V = (fe/Dt) / fx_pix; % scale of the velocity in m/s
% store scales in structure m
m.scale_V = scale_V;
m.ft = fe;
m.fx = fx;
%% Get Histogram displacement 

% Vx = m.Vx(10:end,:,:) - mean(mean(m.Vx(10:end,:,:),2),1);
get_histogram_displacement(m.Vx,W,fig_folder,font_size);

%% Get profile of velocity
i_x = 20;
i_y = 20;
plot_located_profile(m.Vx,i_x,i_y,fx,fe,fig_folder);

disp(mean(abs(m.Vx(i_x,i_y,:))))


%% Get rid off quadratic noise (drone movements)
[nx,ny,nt] = size(m.Vx);
x = (1:1:nx);
y = (ny:-1:1);

% compute the mean component of the velocity field for each frame
Vxmoy = mean(mean(m.Vx,2),1);
Vymoy = mean(mean(m.Vy,2),1);
% reduce the velocity field from its mean value
m.Vx = m.Vx - Vxmoy;
m.Vy = m.Vy - Vymoy;

m.Vx = supress_quadratic_noise(m.Vx,x,y);
m.Vy = supress_quadratic_noise(m.Vy,x,y);


%% Get time Fourier transform
disp('Getting Time Fourier transform')
padding_bool = 1;
add_pow2 = 0;
[FFT_t,TF_spectrum,f] = temporal_FFT(m_scaled.Vx(:,:,:),padding_bool,add_pow2,fe);


%% Try to get wave vector orientation for each frequency

disp('Getting wave vectors')
selected_freq = 0.2; % selected frequencies between which we proceed to the analysis
x_bound = [1 size(m.Vx(:,:,:),1)]; % selected boundaries at which we perform 2D FFT
padding_bool = 1;
add_pow2 = 2; % additional power of 2 for padding 
black_mask = 10;
caxis_amp = -1;

[min_freq, i0] = min(abs(f - selected_freq));
xmin = 1;
xmax = size(FFT_t,1);

disp(i0);
a = FFT_t(xmin:xmax,:,i0); % get the 2D map in temporal Fourier space

% computes the 2D FFT
[shifted_fft,fft_2D,kx,ky] = spatial_FFT(a,padding_bool,add_pow2,fx);
% kx = dimension 1
% ky = dimension 2
% get center position 
center_x = size(shifted_fft,1)/2 + 1;
center_y = size(shifted_fft,2)/2 + 1;

% Mask center 
if black_mask >= 0 % if black_mask < 0, we don't mask the center of the 2D-FFT
    shifted_fft(size(shifted_fft,1)/2+1-black_mask:size(shifted_fft,1)/2+1+black_mask,size(shifted_fft,2)/2+1-black_mask:size(shifted_fft,2)/2+1+black_mask)=zeros;
end

% Try to mask bottom 
%shifted_fft(size(shifted_fft,1)/2 + 30:end,:) = zeros;
figure,
imagesc(kx,ky,abs(shifted_fft)');
%imagesc(abs(shifted_fft));

title(['Frequency : ' num2str(f(i0)) ' Hz'],'Interpreter','latex')
xlabel('$k_x$','Interpreter','latex');
ylabel('$k_y$','Interpreter','latex');
cbar = colorbar();
cbar.Label.String = '$ | \overline{V_x}| (k_x,k_y,f)$';
cbar.Label.Interpreter = 'latex';
if caxis_amp > 0
    caxis([0 caxis_amp])
end 
%axis([-4 4 -4 4])
axis([-1.5 1.5 -1.5 1.5]);
ax = gca;
ax.FontSize = 13;

[Maximum, peak_position] = max(abs(shifted_fft), [],'all', 'linear');
% convert linear index to subscript
sz = size(shifted_fft); % size of the matrix 2D FFT
[row, col] = ind2sub(sz,peak_position); % position of the peak in the matrix
hold on 
% plot position of the peak
kx_peak = kx(row); % dimension 1
ky_peak = ky(col); % dimension 2

plot(kx_peak,ky_peak,'ro');
hold off
%axis([-1.5 1.5 -1.5 1.5]);
%plot(col,row,'ro');

% Scaling of the position
% center_x_scaled and center_y_scaled should be equal to zero
center_x_scaled = kx(center_x);
center_y_scaled = ky(center_y);

%dist(i0-start_idx + 1) = sqrt((row - center_x)^2+(col-center_y)^2);

%distance scaled as inverse of wavelength*2*pi
dist_k = sqrt((kx_peak - center_x_scaled)^2 + (ky_peak - center_y_scaled)^2);
% dist(i0 - start_freq_idx + 1) = dist_k;
disp(['k = ' num2str(dist_k)])
%pause(0.3)

%% Show the demodulated field 
xmin = 1;
xmax = nx;%50*fx;

x = (xmin:1:xmax);
y = (ny:-1:1);
[X,Y] = meshgrid(x,y);

R = real(a(xmin:xmax,:)); % get real intensity map of FFT_temporal at a given frequency

figure,
pcolor(x/fx,y/fx,R')
title(['Frequency : ' num2str(f(i0)) ' Hz'],'Interpreter','latex')
xlabel('$x$ (m)','Interpreter','latex');
ylabel('$y$ (m)','Interpreter','latex');
shading interp
axis([0 size(X,2)/fx 0 size(X,1)/fx])
axis image
cbar = colorbar();
cbar.Label.String = '$ \rm{Re} \left( \overline {V_x}(x,y,f) \right) \: \rm (m.s^{-1})$';
cbar.Label.Interpreter = 'latex';
%     cbar.Label.FontSize = font_size;
if caxis > 0 
    caxis([-caxis_amp caxis_amp])
end 
ax = gca;
ax.FontSize = 13;

%% plot line in the direction of kx and ky

% The line should pass throug this point
x_scaled = x/fx;
y_scaled = y/fx;

[theta,k] = cart2pol(kx_peak,ky_peak);

disp(theta)
x0 = 170;
y0 = 1;
s = (0:-1:-170);

x_line = x0+s*cos(theta);
y_line = y0-s*sin(theta);

figure(11)
hold off
pcolor(x/fx,y/fx,R')
shading interp
hold on 
plot(x_line,y_line,'r--','LineWidth', 3)
axis image
graphe_legende('$x$ (m)','$y$ (m)','$f$ = 0.3 Hz',true)

saveas(gcf,[fig_folder 'champ_demodule_f30cHz_line'],'fig')
saveas(gcf,[fig_folder 'champ_demodule_f30cHz_line'],'png')
saveas(gcf,[fig_folder 'champ_demodule_f30cHz_line'],'eps')


%%
%show velocity field



for i=490:510
    figure(14)
    subplot(1,2,1)
    pcolor(x,y,m.Vx(:,:,i)')
    shading interp
    caxis([-1 1])
    colorbar()
    subplot(1,2,2)
    pcolor(x,y,m.Vy(:,:,i)')
    caxis([-1 1])
    shading interp
    colorbar()
    getframe();
    pause(0.2)
end
%% Interpolate pixels that are closest to this line

% y1 = y(y_line>y_scaled(end));
y1 = round(y_line);
%% Spatial FFT on continuous domain 
x = x_scaled;
y = y_scaled;

[X,Y,T] = ndgrid(x,y,t);

Fx = griddedInterpolant({x,yi,t},m.Vx);
Fy = griddedInterpolant({x,yi,t},m.Vy);


%test 
figure(12);

x0 = 120;
y0 = 36;

[val,i0]= min(abs(x-x0));
[val,j0]= min(abs(yi-y0));
hold off
%plot(Fx(x0*ones(1,nt),y0*ones(1,nt),t))
%hold on
%plot(Fy(x0*ones(1,nt),y0*ones(1,nt),t))

Vx = Fx(x0*ones(1,nt),y0*ones(1,nt),t);
Vy = Fy(x0*ones(1,nt),y0*ones(1,nt),t);

plot(t,-Vx*cos(theta)-Vy*sin(theta))
hold on
plot(t,Vx*sin(theta)-Vy*cos(theta))
hold on
%plot(t,Vx*cos(theta)-Vy*sin(theta))


%%
for k=1:length(s)
    x0 = x_line(k);
    y0 = y_line(k);
    Lx(k,:)=Fx(x0*ones(1,nt),y0*ones(1,nt),t);
    Ly(k,:)=Fy(x0*ones(1,nt),y0*ones(1,nt),t);
    L_ll(k,:)=-Fx(x0*ones(1,nt),y0*ones(1,nt),t)*cos(theta);
end


%plot(Fx(100*ones(1,nt),40*ones(1,nt),t))

figure(13)
pcolor(t,-s,L1);shading interp
graphe_legende('$t$ s)','$L$ (m)','$f$ = 0.3 Hz',true);


%% Do a movie of the velocity field 

save_video = 1;

if save_video
    video_filename = [fig_folder 'Velocity_field_Vy_video.avi']; % folder where the video is saved
    vid = VideoWriter(video_filename);
    vid.FrameRate = 30;
    open(vid)
end 

[nx,ny,nt] = size(m.Vx);
x = (1:1:nx);
y = (ny:-1:1);
[X,Y]=meshgrid(x,y);

velocity_fig = figure;
velocity_fig.Color = [1 , 1 ,1];

for i = 1 : nt
%     V = sqrt(m.Vx(:,:,i).^2 + m.Vy(:,:,i).^2);
    
    pcolor(x/fx,y/fx,m.Vy(:,:,i)')
    shading interp
    title(['Frame ' num2str(i)])
    xlabel('$x$ (m)','Interpreter','latex');
    ylabel('$y$ (m)','Interpreter','latex');
    ax = gca;
    ax.FontSize = 13;
    axis image
    colorbar()
    caxis([-3 3])
    if save_video
        T(i)=getframe(gcf);     
    end

end 

if save_video
    all_valid = true;
    flen = length(T);
    for K = 1 : flen
      if isempty(T(K).cdata)
        all_valid = false;
        fprintf('Empty frame occurred at frame #%d of %d\n', K, flen);
      end
    end
    if ~all_valid
       error('Did not write movie because of empty frames')
    end

    writeVideo(vid,T)
    close(vid)    
end 
