function [lambda,d,freq] = get_attenuation_coef(TF,f,fx,selected_freq,x_bound,freq_thresh,x_bound_thresh,fig_folder,fig_name,save_image,save_video)

% This function computes attenuation coefficient of the wave amplitude along x, for all frequencies
% It returns an array of attenuation coefficient : lambda, and for each
% frequency, the distance of the associated amplitude A_xf to the
% exponential fit. 

% - TF : time Fourier Transform [nx,ny,nt]
% - f : frequency array [nt]
% - fx : scaling for spatial dimension in box/meter
% - selected_freq : 2 x 1 array, with minimum and maximum frequencies we
% study
% - x_bound : min and max indices along x-axis, within which we proceed to
% fitting of the amplitudes
% - freq_thresh : frequency threshold, above which we change the bounds of
% the fitting
% - x_bound_thresh :new min and max boundaries used to fit curves
% associated to frequencies higher than freq_thresh
% - fig_folder : folder where we save plots

min_selected_freq = selected_freq(1);
max_selected_freq = selected_freq(2);

% get indices of frequencies closest to max and min selected frequency
[min_freq, start_freq_idx] = min(abs(f - min_selected_freq)); 
[max_freq, end_freq_idx] = min(abs(f - max_selected_freq));

TF_xf = squeeze(mean(TF(:,:,start_freq_idx:end_freq_idx),2)); % average the TF over y  

% create arrays for plotting
freq = f(start_freq_idx:end_freq_idx);
x_xf = (1:1:size(TF_xf,1))/fx; % in meter

A_xf = abs(TF_xf); % amplitude along x-axis for each frequency
[nx,nf] = size(A_xf);

relevant_indices = start_freq_idx:1:end_freq_idx;

% indices used to fit the exponential decay
i_min = x_bound(1);
i_max = x_bound(2);

i_min_high_freq = x_bound_thresh(1); % lower bound of the spatial domain where the fit is done for high frequencies
i_max_high_freq = x_bound_thresh(2); % upper bound of the spatial domain where the fit is done for high frequencies

lambda = zeros(nf,1); % array of attenuation coefficients
d = zeros(nf,1); % array of distance to plot

decay_fig = figure;
decay_fig.Color = [1,1,1];

if save_video
    video_filename = [fig_folder fig_name '.avi']; % folder where the video is saved
    vid = VideoWriter(video_filename);
    vid.FrameRate = 3;
    open(vid)
end 

for i = 1:nf
%     idx = relevant_indices(i);
    disp(['f = ' num2str(freq(i)) ' Hz'])
    
    if freq(i) > freq_thresh 
        i_min = i_min_high_freq;
        i_max = i_max_high_freq;
    end 
    
    A = A_xf(:,i); % amplitude along x for a given frequency 
    log_A = log10(A); % take the log10 of the amplitude of freq i
    
    % restrict to the boundaries in order to fit
    x_fit = x_xf(i_min:i_max);
    A_fit = log_A(i_min:i_max);
    p = polyfit(x_fit,A_fit,1); % fit log_A by a degree 1 polynome
    lambda(i) = log(10)*p(1); % get attenuation coefficient in m^-1

    disp(['alpha = ' num2str(lambda(i)) ' m'])
    y_poly = 10.^polyval(p,x_fit);
    plot(x_fit,A(i_min:i_max),'o');
    hold on 
    plot(x_fit,y_poly,'r');
    xlabel('$x \: \rm (m)$','Interpreter','latex');
    ylabel('$\langle | V_x | \rangle _y (x,f) \: \rm (m)$','Interpreter','latex');
    if save_video
        title_txt = ['$f = ' num2str(freq(i)) ' \: \rm (Hz)$'];
        title(title_txt,'Interpreter','latex');
    end
    grid on 
    data_txt = 'Data';
    fit_txt = ['$y(x) = C e^{' sprintf('%0.3f',lambda(i)) 'x}$'];
    legend(data_txt,fit_txt,'Interpreter','latex','location','northwest','FontSize',13);
    
    ax = gca;
    ax.FontSize = 13;
    hold off
    
    if save_image
        
        set(gcf,'Units','Inches');
        pos = get(gcf,'Position');
        set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
        print(gcf,'filename','-dpdf','-r0')
        
        frequency_p = strrep(num2str(freq(i)),'.','p'); % replace point by p
        fig_file_fig = [fig_folder 'Attenuation_fit_f_' frequency_p];
        saveas(decay_fig,fig_file_fig,'fig')
    end
    
    if save_video
        T(i)=getframe(gcf);     
    end 
    
    A_red = A(i_min:i_max); % restricted to the region of interest
    d(i) = sum((y_poly - A_red').^2)/sum(A_red.^2); % distance to the fit
    
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
