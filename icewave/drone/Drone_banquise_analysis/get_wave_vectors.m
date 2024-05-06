function [dist,freq] = get_wave_vectors(TF,f,fx,selected_freq,x_bound,padding_bool,add_pow2,black_mask,caxis_amp,fig_folder,vid_name,save_image,save_video)
% This function computes the wave vector associated to a given frequency 
% It also stores a movie of the detected wave vectors 

% It returns : 
% - dist : array of wave vectors, scaled in meter^-1 !
% - freq : array of the studied frequencies

% Function takes as arguments : 
% - TF : time Fourier transform [ny,nx,nt] 
% - f : frequency array 
% - fx : spatial scaling (box / meter)
% - selected_freq : 2 x 1 array, selected frequencies between which we
% compute the wave vectors
% - x_bound : 2 x 1 array, boundaries along the x-axis 
% - padding_bool : boolean to choose whether we pad or not
% - add_pow2 : additional power 2 for padding 
% - caxis_amp : amplitude of the caxis axis 
% - fig_folder : folder where we save plots, images

min_selected_freq = selected_freq(1); % minimal selected frequency on FFT spectrum
max_selected_freq = selected_freq(2); % maximal selected frequency on FFT spectrum

% get indices of frequencies closest to max an min selected frequency
[min_freq, start_freq_idx] = min(abs(f - min_selected_freq)); 
[max_freq, end_freq_idx] = min(abs(f - max_selected_freq));

% get frequencies
freq = f(start_freq_idx:end_freq_idx);

% Values we consider along x-axis to compute the 2D temporal TF field
xmin = x_bound(1);
xmax = x_bound(2);

% Initialize distance array in 2D-Fourier space
dist = zeros(1,start_freq_idx-end_freq_idx);

% Folders where we save images 
new_folder_fig = [fig_folder 'FFT_space/'];
if exist(new_folder_fig,'dir') ~= 7
    mkdir(new_folder_fig)
end

% STARTS

if save_video
    video_filename = [fig_folder vid_name '.avi']; % folder where the video is saved
    vid = VideoWriter(video_filename);
    vid.FrameRate = 3;
    open(vid)
end 
% Loop over all indices
spatial_fourier_fig = figure;
spatial_fourier_fig.Color = [1, 1 ,1];
relevant_indices = start_freq_idx:1:end_freq_idx; 

clear T;
for i = 1:numel(relevant_indices)
    i0 = relevant_indices(i);
    disp(i0);
    a = TF(xmin:xmax,:,i0); % get the 2D map in temporal Fourier space
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
    
    pcolor(kx,ky,abs(shifted_fft)');
    shading interp
    %imagesc(abs(shifted_fft));
    if save_video
        title(['Frequency : ' num2str(f(i0)) ' Hz'],'Interpreter','latex')
    end
    xlabel('$k_x$','Interpreter','latex');
    ylabel('$k_y$','Interpreter','latex');
    cbar = colorbar();
    cbar.Label.String = '$ | \hat{V_x}| (k_x,k_y,f)$';
    cbar.Label.Interpreter = 'latex';
    if caxis_amp > 0
        caxis([0 caxis_amp])
    end 
    %axis([-4 4 -4 4])
    axis([-2.5 2.5 -2.5 2.5]);
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
%     [theta,k_r] = cart2pol(kx_peak,ky_peak); % convert cartesian coordinates to polar
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
%     disp(['k_r = ' num2str(k_r)])
    
    if save_image
        
        set(gcf,'Units','Inches');
        pos = get(gcf,'Position');
        set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
        print(gcf,'filename','-dpdf','-r0')
        
        frequency_p = strrep(num2str(f(i0)),'.','p'); % replace point by p
        fig_file_fig = [new_folder_fig 'FFT_space_f_' frequency_p];
        saveas(spatial_fourier_fig,fig_file_fig,'fig')
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


end 