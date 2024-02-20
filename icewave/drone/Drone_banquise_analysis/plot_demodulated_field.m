function plot_demodulated_field(TF,f,fx,selected_freq,x_bound,caxis_amp,fig_folder,fig_name,save_image,save_video)
% This function plots the real demodulated field of a velocity / height
% field. It also creates a video of the successive demodulated field

% - TF : time-fourier transform [nx,ny,nf] scaled !
% - f : frequency array [nf]
% - selected_freq : 2 x 1 array, frequencies between which we plot the
% demodulated field
% - x_bound : 2 x 1 array, boundaries along x-axis to consider the field to
% be plot


%% Plot a map of the real part of the FFT for each frequency 
% ########### WITH SCALING #############

[nx,ny,nf] = size(TF);
imin = x_bound(1);
imax = x_bound(2);

% create a meshgrid
x = (imin:1:imax);
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

    pcolor(x/fx,y/fx,R')
    shading interp
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