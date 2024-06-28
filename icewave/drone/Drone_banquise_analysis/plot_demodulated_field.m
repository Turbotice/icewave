function plot_demodulated_field(TF,f,fx,selected_freq,x_bound,caxis_amp,left_bool,fig_folder,fig_name,save_image,save_video)
% This function plots the real demodulated field of a velocity / height
% field. It also creates a video of the successive demodulated field

% - TF : time-fourier transform [nx,ny,nf] scaled !
% - f : frequency array [nf]
% - fx : scale in pix/m (or pix/cm depending on which scale you want to
% have along the axis)
% - selected_freq : 2 x 1 array, frequency range for which we plot the
% demodulated field
% - x_bound : 2 x 1 array, boundaries along x-axis to consider the field to
% be plot 
% - caxis_amp : amplitude on the colorbar axis
% - left_bool : boolean to choose how to set direction of x-axis and flip
% image
% - fig_folder : folder where plots and video will be waved
% - fig_name : name under which the video will be saved
% - save_image : boolean to save images or not
% - save_video : boolean to save a video or not


%% Plot a map of the real part of the FFT for each frequency 
% ########### WITH SCALING #############

[nx,ny,nf] = size(TF);
imin = x_bound(1);
imax = x_bound(2);

% create a meshgrid
if left_bool 
    x = (imin:1:imax);
else
    x = (imax:-1:imin);
end 
y = (ny:-1:1);

[X,Y]=meshgrid(x,y);

new_folder = [fig_folder 'FFT_map/'];
if ~exist(new_folder)
    mkdir(new_folder)
end

min_selected_freq = selected_freq(1); % minimal selected frequency on FFT spectrum
max_selected_freq = selected_freq(2); % maximal selected frequency on FFT spectrum

[min_freq, start_freq_idx] = min(abs(f - min_selected_freq)); 
[max_freq, end_freq_idx] = min(abs(f - max_selected_freq));

fig_FFT_map = figure;
fig_FFT_map.Color = [1, 1, 1,];
if save_video
    video_filename = [fig_folder fig_name '.avi']; % folder where the video is saved
    vid = VideoWriter(video_filename);
    vid.FrameRate = 3;
    open(vid)
end 

relevant_indices = start_freq_idx:1:end_freq_idx; 
for i=1:numel(relevant_indices)
    idx = relevant_indices(i);
    disp(idx)
    R = real(TF(imin:imax,:,idx)); % get real intensity map of FFT_temporal at a given frequency
    if left_bool
       R = flip(R,1); 
    end
    pcolor(x/fx,y/fx,R')
    shading interp
%     colormap(redblue)
    xlabel('$x$ (m)','Interpreter','latex');
    ylabel('$y$ (m)','Interpreter','latex');
    shading interp
    axis([0 size(X,2)/fx 0 size(X,1)/fx])
    axis image
    if ~left_bool
       set(gca,'XDir','reverse') 
    end
%     set_Papermode(gcf)
    cbar = colorbar();
    cbar.Label.String = '$ \rm{Re} \left( \hat {V}_x(x,y,f) \right) \: \rm (m.s^{-1})$';
    cbar.Label.Interpreter = 'latex';
%     cbar.Label.FontSize = font_size;
    if caxis_amp > 0 
        caxis([-caxis_amp caxis_amp])
    end 
    ax = gca;
    ax.FontSize = 13;
    
    getframe();
    frequency = f(idx);
    if save_video
        title(['$f = ' num2str(frequency) ' \: \rm(Hz)$'],'Interpreter','latex')
    end 
    set_Papermode(gcf)
    ax = gca;
    ax.FontSize = 13;
    %pause(0.2)
    
    if save_image
        frequency_p = strrep(num2str(frequency),'.','p'); % replace point by p
        fig_file_fig = [new_folder 'FFT_map_f_' frequency_p];
        saveas(fig_FFT_map,fig_file_fig,'fig')
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

end 