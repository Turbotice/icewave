function [S_disp] = extract_branchs_from_Efk(E,f,k,min_prominence,freq_range)

% This function enables to extract points (f,k,A) of main amplitude
% from a space-time spectrum E(f,k)
% Inputs : - E : 2D array (nf,nk), space-time spectrum 
%          - f : 1D array (nf), frequency array 
%          - k : 1D array (nk), wave vector array 
%          - min_prominence, minimal prominence of detected peaks (used on
%          normalize profile)
%          - freq_range : [freq_start freq_end], range of frequencies on
%          which we are detecting branches 

% Outputs : - S_disp : structure containing keys 'f','omega','k','A' of
% each detetected points in the space-time spectrum 

    freq_start = freq_range(1); % start loop at this frequency
    freq_end = freq_range(2); % stop loop at this frequency
    [~,idx_start] = min(abs(freq_start - f));
    [~,idx_end] = min(abs(freq_end - f));
    relevant_idx = (idx_start:1:idx_end);

    figure(12)
    for i = 1 : length(relevant_idx)
        idx = relevant_idx(i);
        disp(['f = ' num2str(f(idx)) ' Hz'])
        max_absolute = max(E(idx,:));
        profile = E(idx,:)/max_absolute; % relative amplitude for a given frequency 

        % detect peaks with a minimum prominence 
        [pks,locs,w,prom] = findpeaks(profile,'MinPeakProminence',min_prominence);
        [y_max,k_peaks] = subpix_precision(profile,k',locs);

        % store detected peaks in a structure
        M_peaks(i).k = k_peaks;
        M_peaks(i).A = y_max*max_absolute;
        M_peaks(i).width = w;
        M_peaks(i).f = f(idx);

        plot(k,profile,'o')
        hold on 
        plot(M_peaks(i).k,y_max,'rd')
        hold off
        xlabel('$k \: \rm (rad.m^{-1})$')
        ylabel('$\hat{V_x}(k,f) \: \rm (m.s^{-1})$')
        pause(0.1)
    end 

    % Creating a structure 
    S_disp = struct('f',[],'k',[],'A',[]);
    for i = 1 : size(M_peaks,2)
        current_f = M_peaks(i).f .* ones(1,length(M_peaks(i).k)); 
        current_k = M_peaks(i).k;
        current_A = M_peaks(i).A;
        S_disp.f = cat(2,S_disp.f,current_f); % store frequency
        S_disp.k = cat(2,S_disp.k,current_k);
        S_disp.A = cat(2,S_disp.A,current_A);
    end 
    S_disp.omega = 2*pi*S_disp.f;  
end 