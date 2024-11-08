function [alpha,C,d] = attenuation_coeff_corrected_direction(FFT_t,f,x,y,facq_x,L0,xmin,f_cutoff,x_cutoff,fig_folder)

% Inputs : 
%   - FFT_t : 3D array, time Fourier transform of velocity field [nx,ny,nf]
%   - f : 1D array, frequencies 
%   - x : 1D array, x-coordinates (in meters)
%   - y : 1D array, y-coordinates (in meters)
%   - facq_x : scalar, acquisition frequency for space coordinates 
%   - L0 : scalar, length over which we look at waves attenuation (in
%   meters)
%   - xmin : scalar, x-coordinate at which we start to fit exponential
%   decay 
%   - f_cutoff : 1D-array, size N, frequencies above which we change range of
%   x-coordinate to perform fitting
%   - x_cutoff : 1D-array, size n, values of xmax associated to each cutoff
%   frequency 
%   - fig_folder : str, name of directory where plots are saved 

% Outputs :
%   - alpha : 1D-array, nf x 1, exponential decaying coefficient of each frequency
%   - C : 1D-array, nf x 1, pre-factor of exponential decay for each frequency
%   - d : 1D-array, nf x 1, norm-2 distance of data points ot fitted
%   exponential curve 

    % parameters used to compute space FFT and binarize space FFT 
    padding_bool = 1; % additional padding for space FFT computation
    add_pow2 = 1; 
    threshold = 0.8; % threshold used to binarize space FFT

    % index of minimal value used to fit data 
    [~,i_min] = min(abs(x - xmin));

    % create empty arrays 
    nf = length(f);
    alpha = zeros(nf,1); % array of attenuation coefficients
    C  = zeros(nf,1); % array of prefactor
    d = zeros(nf,1); % array of distance to plot


    for idx_freq = 1:nf
        current_freq = f(idx_freq);
        freq_txt = replace(num2str(f(idx_freq)),'.','p'); % suffixe used to name figures 
        disp(['f = ' num2str(current_freq) ' Hz'])

        % 2D-FFT of complex field
        field = FFT_t(:,:,idx_freq);
        [shifted_fft,~,kx,ky] = spatial_FFT(field,padding_bool,add_pow2,facq_x);

        % Normalize 2D-FFT
        zmax = max(abs(shifted_fft),[],'all');
        norm_fft = abs(shifted_fft)/zmax;

        peak_idx = 1;
        % binarize the 2D-FFT spectrum
        bin_fft = norm_fft;
        bin_fft(norm_fft > threshold) = 1;
        bin_fft(norm_fft <= threshold) = 0;

        CC = bwconncomp(bin_fft');
        stats = regionprops("table",CC,norm_fft','WeightedCentroid');
        row = round(stats.WeightedCentroid(peak_idx,1));
        col = round(stats.WeightedCentroid(peak_idx,2));
        kx_peak = kx(row);
        ky_peak = ky(col);

        FFT_space_fig = figure(1);
        pcolor(kx,ky,abs(shifted_fft)')
        shading interp
        hold on 
        plot(kx_peak,ky_peak,'ro')
        xlabel('$k_x \: \rm (rad.m^{-1})$')
        ylabel('$k_y \: \rm (rad.m^{-1})$')
        axis([-2 2 -2 2]) 
        hold off

        [theta,~] = cart2pol(kx_peak,ky_peak); % converts to polar coordinates 

        if theta < 0
            % Build several segments for theta < 0
            % initial line 
            x0 = x(1); % minimal value along x axis
            y0 = L0*sin(abs(theta)) + y(1);
            ds = 1/facq_x; % step of curvilinear coordinate
            s = (0:ds:L0); % curvilinear coordinate

            % Points that define the line 
            x_line = x0+s*cos(theta);
            y_line = y0+s*sin(theta);

        else
            % Building several segments for theta > 0
            x0 = (y(end) - L0*sin(theta) - y(1))*tan(theta);
            y0 = y(1);
            ds = 1/facq_x;

            s = (0:ds:L0); % curvilinear coordinate
            x_line = x0 + s*cos(theta);
            y_line = y0 + s*sin(theta);
        end 

        real_field_fig = figure(2); 
        pcolor(x,y,real(field)')
        shading interp
        hold on 
        plot(x_line,y_line,'k--')
        axis image
        xlabel('$x \: \rm (m)$')
        ylabel('$y \: \rm (m)$')
        ax = gca; 
        ax.FontSize = 13;
        set_Papermode(gcf)

        % ask if we choose reoriented field
        prompt = "Do you want to reorient wave field ? y/n [y]: ";
        txt = input(prompt,"s");
        if isempty(txt)
            txt = 'y';
        end

        if txt == 'y'
            % Compute field correctly oriented
            [field_star,x_star,~] = orientation_parallel2propagation(field,theta,1/facq_x,L0);
            hold on 
            plot(x_line,y_line,'k-')
        else
            field_star = field;
            x_star = x;
        end

        hold off
        figname = [fig_folder 'Real_field_f' freq_txt];
        saveas(gcf,figname,'fig')
        saveas(gcf,figname,'pdf')

        % Exponential fit 
        A = squeeze(mean(abs(field_star),2)); % amplitude along x-axis for each frequency

        % indices used to fit the exponential decay

        if isempty(f_cutoff)
            i_max = length(A); % maximum index 
        else 
            idx_cutoff = 0;
            current_fcutoff = f_cutoff(idx_cutoff + 1); % initialize cut off frequency 
            while current_freq > current_fcutoff && idx_cutoff < length(f_cutoff)
                idx_cutoff = idx_cutoff + 1; 
                if idx_cutoff < length(f_cutoff)
                    current_fcutoff = f_cutoff(idx_cutoff + 1);
                end 
            end 

            if idx_cutoff == 0
                i_max = length(A);
            else
                current_xcutoff = x_cutoff(idx_cutoff);
                [~,i_short] = min(abs(x - current_xcutoff));
                i_max = i_short;
            end 
        end 

        decay_fig = figure(3);
        decay_fig.Color = [1,1,1];

        % restrict to the boundaries in order to fit
        x_fit = x_star(i_min:i_max); % in meter !!
        A_red = A(i_min:i_max); % restricted to the region of interest

        log_A = log10(A_red); % take the log10 of the amplitude of freq i
        A_fit = log_A; 
        p = polyfit(x_fit,A_fit,1); % fit log_A by a degree 1 polynome
        alpha(idx_freq) = log(10)*p(1); % get attenuation coefficient in m^-1
        C(idx_freq) = 10.^p(2); % prefactor

        disp(['alpha = ' num2str(alpha(idx_freq)) ' m-1'])
        y_poly = 10.^polyval(p,x_fit);
        plot(x_fit,A_red,'o');
        hold on 
        plot(x_fit,y_poly,'r');
        xlabel('$x \: \rm (m)$','Interpreter','latex');
        ylabel('$\langle | \hat{V_x} | \rangle _y (x,f) \: \rm (m)$','Interpreter','latex');
        grid on 
        data_txt = 'Data';
        fit_txt = ['$y(x) = ' sprintf('%0.2f',C(idx_freq)) ' e^{' sprintf('%0.3f',alpha(idx_freq)) 'x}$'];
        legend(data_txt,fit_txt,'Interpreter','latex','location','northeast','FontSize',13);

        d(idx_freq) = sum((y_poly - A_red').^2)/sum(A_red.^2); % distance to the fit
        disp(['Distance to fit : ' num2str(d(idx_freq)) ''])

        title(['$d = ' sprintf('%0.4f',d(idx_freq)) '$'],'Interpreter','latex')
        ax = gca;
        ax.FontSize = 13;
        set_Papermode(decay_fig);
        hold off

        figname = [fig_folder 'Attenuation' freq_txt];
        saveas(decay_fig,figname,'fig')
        saveas(decay_fig,figname,'pdf')

    end 



end 