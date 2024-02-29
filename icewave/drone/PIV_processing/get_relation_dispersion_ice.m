%% Clear 
clear all;
close all;

%% Loading structure obtained after PIV processing and post-processing
date = '20230310';
% base = '/media/turbots/DATA/thiou/labshared2/Banquise/Rimouski_2023/Data/drone/';
base = 'W:/Banquise/Rimouski_2023/Traitements_donnees/PIV_Sebastien/Data_ice_Dt4_W64/';
% base = 'E:/Stage MIZ/PIVlab_drone';
%base = '/media/turbots/My Passport//Stage MIZ/PIVlab_drone';

% folder = [base '/matData/Total_processed_2023_09_01/'];
folder = base;
prefixe = 'Stephane_PIV_processed_';
suffixe = 'Dt4_b1_DJI_0402_images_post_processed';
filename = [folder prefixe suffixe '.mat'];
%%
disp('Loading Data..');
load(filename);
disp('Data loaded');
%% Creates a folder where to save the generated plots
fig_folder = 'W:/Banquise/Rimouski_2023/Traitements_donnees/PIV_Sebastien/Figures_report/new_scale_ice/';
if ~exist(fig_folder)
    mkdir(fig_folder)
end

%% Scaling 
fe = 30; % Frame rate in Hz
W = 64; % size of the window for DIC algorithm
font_size = 13;

% ##########################################
% Automatize fx_pix !!!!
L_x = 3840; % size of the image in pixel, along x-axis
h_drone = 100.6; % height of the drone in meter
theta_x = 32.75; % semi AFOV of the drone, along x-axis, in °

fx_pix = L_x/(2*h_drone*tan(theta_x*pi/180)); % scale in pixels / meter
fx = fx_pix*2/W;
% ##########################################

% fx = 0.8857; % for W = 64; 
%fx = 0.8857*2; % spatial scale in boxes/meter
Dt = 4; % step between two frames that were compared during the PIV algorithm 

scale_V = (fe/Dt) / fx_pix; % scale of the velocity in m/s
%% Initial picture with scaling 

picture_file = 'W:/Banquise/Rimouski_2023/Data/drone/20230310/contexte/video/DJI_0402_images/im_0001.tiff';
img = imread(picture_file);

%%
% x = (1:1:size(img,2));
% y = (1:1:size(img,1));
new_img = img(:,3000 : end,:);
RI = imref2d(size(new_img));

RI.XWorldLimits = RI.XWorldLimits ./fx_pix; 
RI.YWorldLimits = RI.YWorldLimits ./fx_pix;
figure, imshow(new_img,RI)
xlabel('$x$ (m)','Interpreter','latex');
ylabel('$y$ (m)','Interpreter','latex');
ax = gca;
ax.FontSize = font_size;
% set correctly the image position for a pdf format 
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'filename','-dpdf','-r0')


%%
[ny,nx,nt] = size(m.Vx);

% create a meshgrid
y = (1:1:ny);
x = (nx:-1:1);

[Y,X]=meshgrid(y,x);

% compute the mean component of the velocity field for each frame
% Vxmoy = mean(mean(m.Vx,2),1);
% Vymoy = mean(mean(m.Vy,2),1);
% % reduce the velocity field from its mean value
% m.Vx = m.Vx - Vxmoy;
% m.Vy = m.Vy - Vymoy;

%% Image of the velocity field 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ONLY THE IMAGE IS SCALED %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;

save_boolean = 1;

i = 1;
subplot(2,1,1)
hold off
surf(Y./fx,X./fx,m.Vx(:,:,i)' .*scale_V)
shading interp
view(2)
axis([0 size(X,2)/fx 0 size(X,1)/fx]);
caxis([-1 1].*scale_V)
cbar =colorbar();
cbar.Label.String = '$V_x \: \rm (m.s^{-1})$';
cbar.Label.Interpreter = 'latex';
cbar.FontSize = 13;
%title('$V_x$','Interpreter','latex')
xlabel('$x$ (m)','Interpreter','latex');
ylabel('$y$ (m)','Interpreter','latex');
ax = gca;
ax.FontSize = 13;

subplot(2,1,2)
hold off
surf(Y./fx,X./fx,m.Vy(:,:,i)' .*scale_V)
axis([0 size(X,2)/fx 0 size(X,1)/fx]);
shading interp
view(2)
caxis([-1 1].*scale_V)
cbar =colorbar();
cbar.Label.String = '$V_y \: \rm (m.s^{-1})$';
cbar.Label.Interpreter = 'latex';
cbar.FontSize = 13;
%title('$V_y$','Interpreter','latex')
xlabel('$x$ (m)','Interpreter','latex');
ylabel('$y$ (m)','Interpreter','latex');
ax = gca;
ax.FontSize = 13;

% set correctly the image position for a pdf format 
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'filename','-dpdf','-r0')

fig_file = [fig_folder 'Velocity_map_full_idx_' num2str(i)];
if save_boolean
    saveas(gcf,fig_file,'fig');
    saveas(gcf,fig_file,'pdf');
end


%% Moyenne, ecart type

save = 0;

% computes mean and standard deviation over time 
m.Vxmoy = nanmean(m.Vx,3);
m.Vymoy = nanmean(m.Vy,3);
m.Vxstd = std(m.Vx,[],3);
m.Vystd = std(m.Vy,[],3);
 
fig_average = figure(2);
subplot(2,1,1)
surf(Y./fx,X./fx,m.Vxmoy' .*scale_V)
%title('$\langle V_x \rangle _t$','Interpreter','latex','FontSize',13)
xlabel('$x$ (m)','Interpreter','latex');
ylabel('$y$ (m)','Interpreter','latex');
shading interp
view(2)
axis([0 size(X,2)/fx 0 size(X,1)/fx]);
caxis([-3 3].*scale_V)
cbar = colorbar();
cbar.Label.String = '$\langle V_x \rangle _t \: \rm (m.s^{-1})$';
cbar.Label.Interpreter = 'latex';
cbar.FontSize = 13;
ax = gca;
ax.FontSize = 13;

subplot(2,1,2)
surf(Y./fx,X./fx,m.Vymoy'.*scale_V)
%title('$\langle V_y \rangle _t$','Interpreter','latex')
xlabel('$x$ (m)','Interpreter','latex');
ylabel('$y$ (m)','Interpreter','latex');
shading interp
view(2)
axis([0 size(X,2)/fx 0 size(X,1)/fx]);
caxis([-2 2].*scale_V)
cbar = colorbar();
cbar.Label.String = '$\langle V_y \rangle _t \: \rm (m.s^{-1})$';
cbar.Label.Interpreter = 'latex';
cbar.FontSize = 13;
ax = gca;
ax.FontSize = 13;

% set correctly the image position for a pdf format 
set(fig_average,'Units','Inches');
pos = get(fig_average,'Position');
set(fig_average,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fig_average,'filename','-dpdf','-r0')

fig_std = figure(3);
subplot(2,1,1)
surf(Y./fx,X./fx,m.Vxstd' .*scale_V)
%title('$\sigma (V_x)$','Interpreter','latex')
xlabel('$x$ (m)','Interpreter','latex');
ylabel('$y$ (m)','Interpreter','latex');
shading interp
view(2)
axis([0 size(X,2)/fx 0 size(X,1)/fx]);
caxis([-3 3].*scale_V)
cbar = colorbar();
cbar.Label.String = '$\sigma (V_x) \: \rm (m.s^{-1})$';
cbar.Label.Interpreter = 'latex';
cbar.FontSize = 13;
ax = gca;
ax.FontSize = 13;

subplot(2,1,2)
surf(Y./fx,X./fx,m.Vystd' .*scale_V)
%title('$\sigma (V_y)$','Interpreter','latex')
xlabel('$x$ (m)','Interpreter','latex');
ylabel('$y$ (m)','Interpreter','latex');
shading interp
view(2)
axis([0 size(X,2)/fx 0 size(X,1)/fx]);
caxis([-3 3].*scale_V)
cbar = colorbar();
cbar.Label.String = '$\sigma (V_y) \: \rm (m.s^{-1})$';
cbar.Label.Interpreter = 'latex';
cbar.FontSize = 13;
ax = gca;
ax.FontSize = 13;

% set correctly the image position for a pdf format 
set(fig_std,'Units','Inches');
pos = get(fig_std,'Position');
set(fig_std,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fig_std,'filename','-dpdf','-r0')


fig_file_mean = [fig_folder prefixe suffixe '_velocity_average_full'];
fig_file_std = [fig_folder prefixe suffixe '_velocity_std_full'];

if save 
    saveas(fig_average, fig_file_mean, 'pdf');
    saveas(fig_average, fig_file_mean, 'fig');
    saveas(fig_std, fig_file_std, 'pdf');
    saveas(fig_std, fig_file_std, 'fig');
end


%% Get rid off quadratic noise (drone movements)
[nx,ny,n] = size(m.Vx);
%imax = 105; % maximum index on the x-axis
imax = nx; % maximum idx used for ice 
%imax = nx;
imin = 1; % minimum index on the x-axis
%imin = 25; % minium idx used for waves 
[Y_quad,X_quad] = meshgrid(m.y,m.x(imin:imax));
quadratic_boolean = 1;
% clear Vx_s;
Vx_s = zeros(imax - imin + 1,ny,n);
for i=1:n
    if quadratic_boolean
        Vx = m.Vx(imin:imax,:,i);
        % fit by a quadratic function
        P = polyFit2D(Vx,X_quad,Y_quad,2,2);
        Pth = polyVal2D(P,X_quad,Y_quad,2,2);

        %size(Pth)
        %imagesc(Vx'-Pth')
        Ptot(i,:) = P;
        Vx_s(:,:,i) = Vx-Pth; % get rid off quaratic noise (drone movements)
        %disp('Reduction of quadratic noise');
    else 
        Vx_s(:,:,i) = m.Vx(imin:imax,:,i);
    end
end


%% FFT temporel

% computes the closest power 2 to the signal temporal length
% closest_power_2 = log2(size(m.Vx,3));
% padding_length = floor(power(2,closest_power_2 + 1));
padding_length = 2^nextpow2(size(Vx_s,3));
original_length = size(Vx_s,3);
% computes FFT of Vx
save_boolean = 0;
padding = 1;
% Vx_s = Vx_s - mean(mean(Vx_s,1),2);
Vx_0 = Vx_s - mean(Vx_s,3); % reduction of mean over time
if padding % if we want to pad with zeros
    TF = fft(Vx_0,padding_length,3);
    N = padding_length ;
    disp('Padding used')
else % if we do not pad 
    TF = fft(Vx_0,[],3);
    [nx,ny,nt] = size(TF);
    N = nt;
    disp('No padding')
end

TF = TF/original_length; % Normalization of the fft by the original length of the signal
% get mean vealues of the FFT for each frequency 
TF_inter = squeeze(mean(mean(abs(TF),2),1));

% P2 = abs(TF_t/N);
% P1 = P2(1:N/2+1);
% P1(2:end-1) = 2*P1(2:end-1);

TF_t = TF_inter(1:N/2+1);
TF_t(2:end-1) = 2*TF_t(2:end-1); % multiply by 2 for peaks that are both in positive an negative frequencies

%fe = 30;%to be modified
f = fe*(0:(N/2))/N;
% f = linspace(0,fe/2,N/2+1);
%%
% plot FFT spectrum, fe needs to be correctly chosen 
fig_FFT_scaled = figure(6);
loglog(f,TF_t .*scale_V)
xlabel('$f$ (Hz)','Interpreter','latex')
ylabel('$\overline {\langle V_x \rangle _{x,y}} (f) \: \rm (m.s^{-1})$','Interpreter','latex')
ax = gca;
ax.FontSize = font_size;
grid on
axis([0.001 50 1/1e6 0.01])

% set correctly the image position for a pdf format 
set(fig_FFT_scaled,'Units','Inches');
pos = get(fig_FFT_scaled,'Position');
set(fig_FFT_scaled,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fig_FFT_scaled,'filename','-dpdf','-r0')

% % plot FFT array 
% fig_FFT = figure(7);
% loglog(TF_t)
% xlabel('FFT index')
% ylabel('Amplitude')
% title('Temporal FFT')
% grid on
% hold on
%fig_file = [fig_folder prefixe suffixe '_FFT'];
fig_file_2 = [fig_folder prefixe suffixe '_FFT_scaled_ice_quadratic'];
% saving graphs
if save_boolean 
%     saveas(fig_FFT, fig_file, 'fig');
%     saveas(fig_FFT, fig_file, 'pdf');
    saveas(fig_FFT_scaled, fig_file_2, 'fig');
    saveas(fig_FFT_scaled, fig_file_2, 'pdf');
end 

%%

figure;
loglog(f,TF_t_quad .*scale_V);
hold on 
loglog(f,TF_t .*scale_V);
xlabel('$f$ (Hz)','Interpreter','latex')
ylabel('$\overline {\langle V_x \rangle _{x,y}} (f) \: \rm (m.s^{-1})$','Interpreter','latex')
legend('quadratic corrections','no corrections','Interpreter','latex','FontSize', font_size,'location','southwest')
grid on
ax = gca;
ax.FontSize = font_size;
axis([0.001 50 1/1e6 0.02])

set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'filename','-dpdf','-r0')

%% Plot a map of the real part of the FFT for each frequency 
% ########### WITH SCALING #############

[ny,nx,nf] = size(TF);
imax = ny; % maximum index on the x-axis
%imax = 105;
imin = 1;
% create a meshgrid
y = (imin:1:imax);
x = (nx:-1:1);

[Y,X]=meshgrid(y,x);

new_folder = [fig_folder 'FFT_map_ice/'];
if ~exist(new_folder)
    mkdir(new_folder)
end

new_folder_fig = [new_folder 'fig/'];
if ~exist(new_folder_fig)
    mkdir(new_folder_fig)
end
new_folder_pdf = [new_folder 'pdf/'];
if ~exist(new_folder_pdf)
    mkdir(new_folder_pdf)
end

% parameters for saving 
save_image = 0;
save_gif = 0;
save_video = 1;

min_selected_freq = 0.05; % minimal selected frequency on FFT spectrum
max_selected_freq = 0.7; % maximal selected frequency on FFT spectrum

[min_freq, start_freq_idx] = min(abs(f - min_selected_freq)); 
[max_freq, end_freq_idx] = min(abs(f - max_selected_freq));

gif_filename = [fig_folder 'real_FFT_map.gif'];

fig_FFT_map = figure;
if save_video
    video_filename = [fig_folder prefixe 'Reconstruction_FFT_temporal_ice.avi']; % folder where the video is saved
    vid = VideoWriter(video_filename);
    vid.FrameRate = 3;
    open(vid)
end 
relevant_indices = start_freq_idx:1:end_freq_idx; 
for i=1:numel(relevant_indices)
    idx = relevant_indices(i);
    disp(idx)
    R = real(TF(:,:,idx)); % get real intensity map of FFT_temporal at a given frequency

    surf(Y/fx,X/fx,R' .*scale_V)
    xlabel('$x$ (m)','Interpreter','latex','FontSize', font_size);
    ylabel('$y$ (m)','Interpreter','latex','FontSize', font_size);
    shading interp
    view(2)
    axis([0 size(X,2)/fx 0 size(X,1)/fx])
    axis image
    cbar = colorbar();
    cbar.Label.String = '$ \rm{Re} \left( \overline {V_x}(x,y,f) \right) \: \rm (m.s^{-1})$';
    cbar.Label.Interpreter = 'latex';
%     cbar.Label.FontSize = font_size;
    caxis([-150 150] .*scale_V /original_length)
    ax = gca;
    ax.FontSize = font_size;
    
    getframe();
    frequency = f(idx);
    if save_video
        title(['Frequency : ' num2str(frequency) ' (Hz)'],'Interpreter','latex')
    end 
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(gcf,'filename','-dpdf','-r0')
    %pause(0.2)
    
    if save_image
        frequency_p = strrep(num2str(frequency),'.','p'); % replace point by p
        fig_file_fig = [new_folder_fig 'FFT_map_f_' frequency_p];
        fig_file_pdf = [new_folder_pdf 'FFT_map_f_' frequency_p];
        saveas(fig_FFT_map,fig_file_fig,'fig')
        saveas(fig_FFT_map,fig_file_pdf,'pdf')
    end
    
    % save a gif 
    if save_gif
        
        frame = getframe(fig_FFT_map);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);

        if i == 1 %starting number of loopp
            imwrite(imind,cm,gif_filename,'gif','LoopCount',inf);
        else
            imwrite(imind,cm,gif_filename,'gif','WriteMode','append');
        end
    end 
    
    if save_video
        T(i)=getframe(gcf);     
    end 
    
end

% Write video
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

%% Do fft-2D on ALL frequencies 
% ##############################

min_selected_freq = 0.05; % minimal selected frequency on FFT spectrum
max_selected_freq = 0.65; % maximal selected frequency on FFT spectrum

% get indices of frequencies closest to max an min selected frequency
[min_freq, start_freq_idx] = min(abs(f - min_selected_freq)); 
[max_freq, end_freq_idx] = min(abs(f - max_selected_freq));

% get frequencies
freq = f(start_freq_idx:end_freq_idx);
omega = 2*pi*freq;

% Values we consider on the 2D temporal TF field
xmin = 1;
% xmin = 25;
% xmax = 105;
xmax = size(TF,1);

% Initialize distance array in 2D-Fourier space
dist = zeros(1,start_freq_idx-end_freq_idx);

% parameters for saving 
save_image = 0;
save_gif = 0;
save_video = 1;

% Folders where we save images 
new_folder_fig = [fig_folder 'FFT_space/fig/'];
if ~exist(new_folder_fig)
    mkdir(new_folder_fig)
end
new_folder_pdf = [fig_folder 'FFT_space/pdf/'];
if ~exist(new_folder_pdf)
    mkdir(new_folder_pdf)
end

% STARTS

if save_video
    video_filename = [fig_folder prefixe 'Spatial_Fourier_space_ice.avi']; % folder where the video is saved
    vid = VideoWriter(video_filename);
    vid.FrameRate = 3;
    open(vid)
end 
% Loop over all indices
spatial_fourier_fig = figure(32);
relevant_indices = start_freq_idx:1:end_freq_idx; 

clear T;
for i = 1:numel(relevant_indices)
    i0 = relevant_indices(i);
    disp(i0);
    a = TF(xmin:xmax,:,i0); % get the 2D map in temporal Fourier space
    %a = a - mean(mean(a,2),1);
    
    [Nx, Ny] = size(a);
    Nb_elements = Nx*Ny; % number of elements of matrix a
    % padding
    add_pow2 = 3; % additional power 2  for padding
    padding_x = 2^(nextpow2(Nx) + add_pow2);
    padding_y = 2^(nextpow2(Ny) + add_pow2);
    % get center position 
    center_x = padding_x/2 + 1;
    center_y = padding_y/2 + 1;
    % computes the 2D FFT
    fft_2D = fft2(a,padding_x,padding_y);
    fft_2D = fft_2D/Nb_elements; % normalization by the initial number of elements
    % #################################################
    % SCALING : arrays used for plot of the spacial FFT
    % #################################################
    kx = 2*pi*fx*(-padding_x/2:padding_x/2-1)/padding_x;
    ky = 2*pi*fx*(-padding_y/2:padding_y/2-1)/padding_y;
    
    % find the maximum
    shifted_fft = fftshift(abs(fft_2D));
    black_mask = 10; % mask used to get rid off the center for the ice
    shifted_fft(size(shifted_fft,1)/2+1-black_mask:size(shifted_fft,1)/2+1+black_mask,size(shifted_fft,2)/2+1-black_mask:size(shifted_fft,2)/2+1+black_mask)=zeros;
    % Try to mask bottom 
    %shifted_fft(size(shifted_fft,1)/2 + 30:end,:) = zeros;
    
    imagesc(kx,ky,abs(shifted_fft)' .*scale_V);
    %imagesc(abs(shifted_fft));
    if save_video
        title(['Frequency : ' num2str(f(i0)) ' Hz'],'Interpreter','latex')
    end
    xlabel('$k_x$','Interpreter','latex');
    ylabel('$k_y$','Interpreter','latex');
    [Maximum, peak_position] = max(shifted_fft, [],'all', 'linear');
    cbar = colorbar();
    cbar.Label.String = '$ | \overline{V_x}| (k_x,k_y,f)$';
    cbar.Label.Interpreter = 'latex';
    caxis([0 0.0005])
    %axis([-4 4 -4 4])
    axis([-1.5 1.5 -1.5 1.5]);
    ax = gca;
    ax.FontSize = 13;
    % convert linear index to subscript
    sz = [padding_x padding_y]; % size of the matrix 2D FFT
    [row, col] = ind2sub(sz,peak_position); % position of the peak in the matrix
    hold on 
    % plot position of the peak
    kx_peak = kx(row);
    ky_peak = ky(col);
    
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
    dist(i0 - start_freq_idx + 1) = dist_k;
    disp(['k = ' num2str(dist_k)])
    %pause(0.3)
    
    if save_image
        
        set(gcf,'Units','Inches');
        pos = get(gcf,'Position');
        set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
        print(gcf,'filename','-dpdf','-r0')
        
        frequency_p = strrep(num2str(f(i0)),'.','p'); % replace point by p
        fig_file_fig = [new_folder_fig 'FFT_space_f_' frequency_p];
        fig_file_pdf = [new_folder_pdf 'FFT_space_f_' frequency_p];
        saveas(spatial_fourier_fig,fig_file_fig,'fig')
        saveas(spatial_fourier_fig,fig_file_pdf,'pdf')
    end
    
    % save a gif 
    if save_gif
        
        frame = getframe(spatial_fourier_fig);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);

        if i == 1 %starting number of loopp
            imwrite(imind,cm,gif_filename,'gif','LoopCount',inf);
        else
            imwrite(imind,cm,gif_filename,'gif','WriteMode','append');
        end
    end 
    % create a video
    if save_video
        T(i)=getframe(gcf);     
    end 
    
end

% Write video
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

%% Plot dispersion relation using 2D-FFT

new_f = freq(1:end);
new_omega = omega(1:end);
wave_vector = dist;

% first plot before masking
figure(15)
loglog(wave_vector,omega,'o');
grid on 
%% Computing filtered datas
% masking some datas
%mask = (wave_vector > 0.12) & (omega < 2.78);
opposite_mask = (omega > 2.78) & (wave_vector < 0.13);
mask = ~opposite_mask;
fitted_omega = new_omega(mask);
fitted_k = wave_vector(mask);
%%
opposite2 = (fitted_omega < 0.65) & (fitted_k > 0.16) ;
mask2 = ~opposite2;
%mask2 = (fitted_omega > 0.59) ;

fitted_omega = fitted_omega(mask2);
fitted_k = fitted_k(mask2);
%% 
opposite3 = (fitted_k < 0.12) & (fitted_omega > 0.65);
mask3 = ~opposite3;
fitted_omega = fitted_omega(mask3);
fitted_k = fitted_k(mask3);
%%
mask4 = fitted_omega > 0.4;
fitted_omega = fitted_omega(mask4);
fitted_k = fitted_k(mask4);
%Saving relevant data
save_boolean = 1;

if save_boolean
    save_file_name = [fig_folder 'Relevant_datas_waves_2023_01_04_Elie_scaling_addpad3.mat'];
    clear save
    disp('Loading plot data');
    save(save_file_name,'fitted_omega','fitted_k'); 
    disp('Loading plot data, DONE.');
end
%% Final Plot !

data_plot = [fig_folder 'Relevant_datas_ice_2023_12_04_Elie_scaling_addpad3.mat'];
load_bool = 1;
if load_bool
    load(data_plot)
    disp('Data to plot loaded');
end 

save_boolean = 0;

%Physical_parameters
g = 9.81; % intensity of gravity
h = 2.0; % water depth in meter
D = 3.2*10^7; % Flexion modulus
rho = 1000; % water density

% minimal values of the plot on x and y axis
xlim_min = 0.1;
xlim_max = 8;
ylim_min = 0.5;
ylim_max = 10;

% Computation of the different models 
k_list = linspace(xlim_min,xlim_max,100); % array of k for plots
deep_water = sqrt(g*k_list);
shallow_water = sqrt(g*h*k_list.^2);
flexural = sqrt(D*(k_list.^5)/rho);
gravito_waves = sqrt(g*k_list .* tanh(k_list*h));

% Power law fitting 
% data to plot
mask_powerlaw = (fitted_omega < 2.0) & (fitted_k > 0.53) ;
k_powerlaw = fitted_k(mask_powerlaw);
omega_powerlaw = fitted_omega(mask_powerlaw);
l1 = fminsearch(@(s)powerfit(k_powerlaw,omega_powerlaw,s),[1,1]);
% l1 = fminsearch(@(s)powerfit(fitted_k,fitted_omega,s),[1,1]);
yth = powerfun(k_list,l1); % fitted exponential function

% Plot
relation_disp_scaled = figure(34);
%loglog(wave_vector,new_omega,'o');
loglog(fitted_k,fitted_omega,'o');
hold on 
loglog(k_list,deep_water,'k-');
hold on 
loglog(k_list,shallow_water,'b-');
% hold on 
% loglog(k_list,yth,'r-');
hold on 
loglog(k_list,gravito_waves,'k--');
% hold on 
% loglog(k_powerlaw,omega_powerlaw,'o','MarkerFaceColor',[0.4660 0.6740 0.1880],'MarkerEdgeColor','black');
title('Dispersion Relation');
xlabel('Wave number $k$ $(m^{-1})$','Interpreter','latex');
ylabel('Pulsation $\omega$ $(s^{-1})$','Interpreter','latex');
axis([xlim_min xlim_max ylim_min ylim_max]);
shallow_txt = ['Shallow, $h = ' sprintf('%0.1f',h) '\: \rm m$'];
power_law_txt = ['$\omega = ' sprintf('%0.2f',l1(2)) 'k^{' sprintf('%0.2f',l1(1)) '}$'];
gravito_txt = ['$\omega = \sqrt{(g  k  tanh(k h))}$'];
legend('Data','Deep water', shallow_txt, gravito_txt, 'Interpreter','latex','Location','northwest');
grid on

% Saving the plot
%fig_folder = 'C:/Users/sebas/Stage_MIZ/Traitement_PIV_Rimouski_20230310/';
file_relation_disp_scaled = [fig_folder 'Ice_Relation_disp_scaled_Elie_scaling_addpad3'];
if save_boolean 
    saveas(relation_disp_scaled,file_relation_disp_scaled,'pdf');
    saveas(relation_disp_scaled,file_relation_disp_scaled,'fig');
end

%% Try to get the attenuation coefficient 

min_selected_freq = 0.05; % minimal selected frequency on FFT spectrum
max_selected_freq = 2.0; % maximal selected frequency on FFT spectrum

% get indices of frequencies closest to max and min selected frequency
[min_freq, start_freq_idx] = min(abs(f - min_selected_freq)); 
[max_freq, end_freq_idx] = min(abs(f - max_selected_freq));

TF_xf = squeeze(mean(TF(:,:,start_freq_idx:end_freq_idx),2)); % average the TF over y  

% create arrays for plotting
freq_xf = f(start_freq_idx:end_freq_idx);
x_xf = (1:1:size(TF_xf,1));

[X_xf,F_xf] = meshgrid(x_xf,freq_xf); 

A_xf = abs(TF_xf);
figure, surf(X_xf,F_xf,A_xf')
shading interp
view(2)
xlabel('$x \: \rm (m)$','Interpreter','latex');
ylabel('$f \: \rm (Hz)$','Interpreter','latex');
ax = gca;
ax.FontSize = 13;

%% Fit each amplitude A_x of each frequency 

min_selected_freq = 0.1; % minimal selected frequency on FFT spectrum
max_selected_freq = 1.5; % maximal selected frequency on FFT spectrum

% get indices of frequencies closest to max an min selected frequency
[min_freq, start_freq_idx] = min(abs(freq_xf - min_selected_freq)); 
[max_freq, end_freq_idx] = min(abs(freq_xf - max_selected_freq));
[nx,nf] = size(A_xf(:,start_freq_idx:end_freq_idx));

freq = freq_xf(start_freq_idx:end_freq_idx);
relevant_indices = start_freq_idx:1:end_freq_idx;
% indices used to fit the exponential decay
i_min = 60;
i_max = nx;

freq_thresh = 0.46; % frequency at which we change the domain where the power law is fitted
i_min_high_freq = 90; % lower bound of the spatial domain where the fit is done for high frequencies

lambda = zeros(nf,1);
d = zeros(nf,1);
decay_fig = figure;
for i = 1:nf
    idx = relevant_indices(i);
    
    if freq(i) > freq_thresh 
        i_min = i_min_high_freq;
    end 
    
    A = A_xf(:,idx);
    log_A = log10(A); % take the log10 of the amplitude of freq i

        x_fit = x_xf(i_min:i_max);
        A_fit = log_A(i_min:i_max);
    p = polyfit(x_fit,A_fit,1); % fit log_A by a degree 1 polynome
    lambda(i) = p(1); % get attenuation coefficient 

    y_poly = 10.^polyval(p,x_fit);
    plot(x_fit,A(i_min:i_max));
    hold on 
    plot(x_fit,y_poly,'r');
    xlabel('$x \: \rm (m)$','Interpreter','latex');
    ylabel('$\langle | V_x | \rangle _y (x,f) \: \rm (m)$','Interpreter','latex');
    title_txt = ['$f = ' num2str(freq(i)) ' \: \rm (Hz)$'];
    title(title_txt,'Interpreter','latex');
    ax = gca;
    ax.FontSize = 13;
    hold off
    
    
    A_red = A(i_min:i_max);
    d(i) = sum((y_poly - A_red').^2)/sum(A_red.^2);
    pause(0.3)
end 
    

%% Plot evolution of the attenuation coefficient with the frequency 

% masking to compute a first power law with well fitted datas
mask_law = d'<0.05; 
freq_mask = freq(mask_law);
lambda_mask = lambda(mask_law);

mask_fit = freq_mask < 0.6;
freq_fit = freq_mask(mask_fit);
lambda_fit = lambda_mask(mask_fit);
l1 = fminsearch(@(s)powerfit(freq_fit,lambda_fit',s),[1,1]);
% l1 = fminsearch(@(s)powerfit(fitted_k,fitted_omega,s),[1,1]);
f_list = linspace(0.1,2,100);
yth = powerfun(f_list,l1); % fitted exponential function

% Take more datas and get rid off the ones too far from the first fitted power
% law 

mask_secondary = d'<0.2; % secondary mask used to plot more datas
freq_secondary = freq(mask_secondary);
lambda_secondary = lambda(mask_secondary);
theory = powerfun(freq_secondary,l1);
dist_to_fit = sum((theory - lambda_secondary).^2)/sum(lambda_secondary.^2);

mask = dist_to_fit < 1;
freq_plot = freq_secondary(mask);
lambda_plot = lambda_secondary(mask);

attenuation_fig = figure;
loglog(freq_secondary,lambda_secondary,'o');
hold on 
plot(f_list,yth,'r--');
hold on 
loglog(freq_plot,lambda_plot,'dr');
% plot(freq_secondary,theory,'k-')
grid on 
xlabel('$f \: \rm (Hz)$','Interpreter','latex');
ylabel('$\alpha \: \rm (m^{-1})$','Interpreter','latex');
axis([0.1 2 0.0001 1])
ax = gca;
ax.FontSize = 13;

power_law_txt = ['$\alpha = ' sprintf('%0.2f',l1(2)) 'f^{' sprintf('%0.2f',l1(1)) '}$'];
legend('Data',power_law_txt,'Interpreter','latex','Location','southeast')
%% FFT temporal smooth 

save_boolean = 1;
% create cuts in the matrix
[nx,ny,nt] = size(Vx_s);
nb_cuts = 60;
r = rem(nt,nb_cuts); % reminder
q = (nt - r)/nb_cuts; % quotient 

clear Vx_cut
clear TF_cut_t
%if r < q/2 % last cut is too small to be considered
    for j = 1 : nb_cuts 
        Vx_cut(:,:,:,j) = Vx_s(:,:,(j-1)*q + 1 : j*q);
    end
    
% else % last cut is sufficiently big to be considered
%     nb_cuts = nb_cuts + 1;
%     for j = 1 : nb_cuts
%         Vx_cut(:,:,:,j) = Vx_s(:,:,(j-1)*q + 1 : j*q);
%     end 

% perform FFT temporal 
for j = 1:nb_cuts
    disp(j);
    signal_length = size(Vx_cut(:,:,:,j),3);
    addpad = 3;
    padding_length = 2^(nextpow2(signal_length)  + addpad);

    TF_cut_t(:,:,:,j) = fft(Vx_cut(:,:,:,j),padding_length,3)/signal_length;
end 

TF_moy = squeeze(mean(mean(mean(abs(TF_cut_t),4),2),1)); % average over all cuts and over space
N = padding_length;

TF_smooth = TF_moy(1:N/2+1);
TF_smooth(2:end-1) = 2*TF_smooth(2:end-1); % multiply by 2 for peaks that are both in positive an negative frequencies

f = fe*(0:(N/2))/N;

% plot FFT spectrum, fe needs to be correctly chosen 
fig_FFT_scaled = figure(6);
loglog(f,TF_smooth)
xlabel('Frequency (Hz)')
ylabel('Amplitude')
title('Temporal FFT')
grid on
hold on 

% plot FFT array 
% fig_FFT = figure(7);
% loglog(TF_t)
% xlabel('FFT index')
% ylabel('Amplitude')
% title('Temporal FFT')
% grid on
% hold on
%fig_file = [fig_folder prefixe suffixe '_FFT'];
fig_file_2 = [fig_folder prefixe suffixe '_FFT_smooth_scaled'];
% saving graphs
if save_boolean 
%     saveas(fig_FFT, fig_file, 'fig');
%     saveas(fig_FFT, fig_file, 'pdf');
    saveas(fig_FFT_scaled, fig_file_2, 'fig');
    saveas(fig_FFT_scaled, fig_file_2, 'pdf');
end 

%% Definition of the useful functions to fit the amplitude of each Fourier componant along the horizontal

    function yth=expfun(x,l)
        yth=l(1)*exp(-l(2)*x);
    end

    function d = expfit(x,y,l)
        yth = expfun(x,l);
        d = sum((yth-y).^2);
    end

    function yth=cosfun(x,l)
        yth=sin(l(1)*x+l(2));
    end

    function yth=sinexpfun(x,l)
       yth =  expfun(x,[l(1) l(2)]).*cosfun(x,[l(3) l(4)]);
    end

    function d = sinexpfit(x,y,l)
        yth = sinexpfun(x,l);
        d = sum((yth-y).^2);
    end

    function d = sinfit(x,y,l1,l2)
        yth = expfun(x,l1).*cosfun(x,l2);
        d = sum((yth-y).^2);
    end
    
    function yth = powerfun(x,l)
        yth = l(2)*x.^l(1);
    end 

    function d = powerfit(x,y,l)
        yth = powerfun(x,l);
        d = sum((yth-y).^2);
    end 