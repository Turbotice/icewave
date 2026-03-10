
%% inputs of the program

%date = '20241105';
date = '0506';
acq_num = '2';
camera_SN = '40300722';%'40437120';%
%daily_folder = ['/media/turbots/DATA/thiou/storageshared/Banquise/Vasco/Frigo_pmmh/' date '/manip_f_k_membrane/Acquisition_1/camera_40300722/'];% pour l'instant suppose qu'il n'y a qu'une manip par jour
%daily_folder = ['F:/manip_grenoble2024/manips_relation_dispersion/' date '/Acquisition_' acq_num '/camera_' camera_SN '/'];
%daily_folder = ['G:/Grenoble/' date '/manip_relation_dispersion/Acquisition_' acq_num '/camera_' camera_SN '/'];
daily_folder = ['R:/Gre25/' date '/cameras/manip_relation_dispersion/Acquisition_' acq_num '/camera_' camera_SN '/'];

% Get a list of all files and folders in the directory
items = dir([daily_folder '*Hz']);

% Initialize an empty cell array to hold folder names
folderList = {};

% Loop through the items to identify folders
for i = 1:length(items)
    % Check if it's a folder and not '.' or '..'
    if items(i).isdir && ~strcmp(items(i).name, '.') && ~strcmp(items(i).name, '..')
        folderList{end+1} = items(i).name; % Add folder name to the list
    end
end
folderListstr = string(folderList);

% Initialize two empty arrays for the two lists of floats
list_f_exc_str = {};
list_freq_acq_str = {};
%list_f_exc = [];
%list_freq_acq = [];

% Loop through each string in the input list
for i = 1:length(folderListstr)
    % Split the string using the underscore as delimiter
    splitStr = strsplit(folderListstr{i}, '_');
    
    % Extract the first and second parts, remove the 'Hz' suffix, and convert to numbers
    f_exc_temp = erase(splitStr{1}, 'Hz');
    freq_acq_temp = erase(splitStr{2}, 'Hz');
    
    % Append the numbers to the respective lists
    list_f_exc_str{end+1} = f_exc_temp;
    list_freq_acq_str{end+1} = freq_acq_temp;
end



%% Definition of the folder to process and where to save the results
for index_freq = 1:length(list_f_exc_str)
    f_exc_str = list_f_exc_str{index_freq};
    freq_acq_str = list_freq_acq_str{index_freq};
    
    f_exc = str2double(f_exc_str);
    freq_acq = str2double(freq_acq_str);

    optional_intermediate_dir = '';
    base = [daily_folder '/' optional_intermediate_dir '/' f_exc_str 'Hz_' freq_acq_str 'Hz/'];
    %base = ['/media/turbots/DATA/thiou/storageshared/Banquise/Vasco/Frigo_pmmh/' date '/' num2str(f_exc,8) 'Hz_' num2str(freq_acq,8) 'Hz/'];
    
    %base = ['/media/turbots/DATA/thiou/storageshared/Banquise/Vasco/Frigo_pmmh/' date '/test_sans_glace/' num2str(f_exc,8) 'Hz_' num2str(freq_acq,8) 'Hz/'];
    
    %base = ['X:/Banquise/Vasco/Frigo_pmmh/' date '/' num2str(f_exc,8) 'Hz_' num2str(freq_acq,8) 'Hz/'];
    % base = ['X:/Banquise/Vasco/Frigo_pmmh/04042024/' num2str(f_exc,8) 'Hz_' num2str(freq_acq,8) 'Hz/'];
    % base = replace(base,'.','p');
    
    %base = ['/media/turbots/DATA/thiou/storageshared/Banquise/Vasco/Frigo_pmmh/' date '/' optional_intermediate_dir '/sweep/'];
    
    % folder = [base date '/contexte/video/'];
    optional_sufix = '';
    folder = [base 'images' optional_sufix '/']; 
    filelist = dir([folder 'Basler*.tiff']); dirnames={};
    [dirnames{1:length(filelist),1}] = deal(filelist.name);
    dirnames = sortrows(dirnames);
    amount = length(dirnames);
    disp(amount)
    %save_folder = '/media/turbots/DATA/thiou/labshared2/Banquise/Rimouski_2023/Traitement_donnees/PIV_Sebastien';
    save_folder = [base 'matData' optional_sufix '/'];
    %dirsave=[save_folder 'raw_datas/'];
    dirsave = save_folder;
    
    for i=1:length(filelist)
        disp([num2str(i) ' : ' filelist(i).name])
    end
    
    %% Main script
    
    if ~exist(dirsave)
        mkdir(dirsave)
    end
    %
    
    % pause(10000)
    for i=1:1%amount
        folderfile = filelist(i).folder;
        disp('checkpoint 1')
        disp(folderfile)
        name=filelist(i).name;%strcat(folder, '/', dirnames{i,1},'/cam_2/');
        disp('checkpoint 2')
        disp(name)
        %    jetstring=extractBetween(dirnames{i}, 'jet_', '_f');
        %    jet=strcmp(extractBetween(dirnames{i}, 'jet_', '_f'), 'on');
        %    freq=extractBetween(name, strcat(jetstring, '_f'), 'Hz')
        %    f=str2double(strrep(freq '_', '.'))
        
        
        forcing='';%extractBetween(dirnames{i}, 'mikro50mm_', '_z');
    %     namesave= extractBetween([name '.'],'_offset','.');
        namesave = dirnames{i,1};
        disp('checkpoint 3')
        %namesave = num2str(i);
        %disp(forcing)
        
        %    matname=string(strcat(['D:\Surface waves\Wave jet interaction\' date '\matData\'],'jet_', jetstring,'_f_', freq,'Hz_down.mat'));
            %namesave=extractBetween(name, 'z740mm_', '_1');
            
        i0 = 0; %process starting from image i0
        N = 0; %number of frames to analyze
        Dt = 50; %ratio between the fps and the scanning frequency (number of image between image A and image B)
        b = 1; %number of images between image A and image A' (from one step to an other)
        W = 64;
        
        %xROI = 1 ;
        %widthroi = 3840;
        %yROI = 1;
        %heightroi = 2160;
    
        %namesave = [namesave 'N_' num2str(N) '_i0_' num2str(i0)];
        prefix = ['PIV_processed_i0' num2str(i0) '_N' num2str(N) '_Dt' num2str(Dt) '_b' num2str(b) '_W' num2str(W) '_full'];
        matname=string(strcat(dirsave,prefix,'.mat'));
        disp(matname)
        if ~isfile(matname)
            disp('checkpoint 4')
            disp(name)
            fullname = [folderfile '/']; % '/' name];
            disp('checkpoint 5')
            disp(fullname)
    %        [u,v,s,p] = PIVlab_commandline_movie_banquise(fullname,'',N,i0,Dt,b);
            [u,v,s,p,xpix,ypix]=PIVlab_commandline_parallel_stereo(fullname,'',N,i0,Dt,b,W);
            disp('checkpoint 6')
            %[u,v,s,p]=PIVlab_commandline(name);
            
            % save variables u , v ,s and p in a .mat file
            save(matname,'u','v','s','p')
            
        else
            disp(matname)
            disp('Mat file already exist, change matname to run')
        end
        disp(i)
    end
    %% debut code suivant pour post traitement de la reconstruction
    %%
    %##########################################################################
    % This code is used to postprocess a data file from a PIV analysis computed
    % thanks to PIVlab. It also enables us to generate a structure from the
    % postprocessed datas and save it under a .mat file
    %##########################################################################
    
    %% Post-process raw datas from PIV analysis and create an associated structure 
    %clear all
    %base = ['/media/turbots/DATA/thiou/storageshared/Banquise/Vasco/Frigo_pmmh/29032024/200Hz/matData/'];
    %base = ['/media/turbots/DATA/thiou/storageshared/Banquise/Vasco/Frigo_pmmh/' date '/test_sans_glace/' num2str(f_exc,8) 'Hz_' num2str(freq_acq,8) 'Hz/matData/'];
    %base = ['/media/turbots/DATA/thiou/storageshared/Banquise/Vasco/Frigo_pmmh/' date '/' num2str(f_exc,8) 'Hz_' num2str(freq_acq,8) 'Hz/matData/'];
    
    base = [base 'matData' optional_sufix '/'];
    %base = ['/media/turbots/DATA/thiou/storageshared/Banquise/Vasco/Frigo_pmmh/02042024/tests_avec_membrane/50Hz_b/matData/'];
    folder = base; %'matData/waves_003/'];% folder of raw datas
    filename = prefix;%'PIV_processed_i00_N0_Dt1_b1_W16_full.mat';
    fullname = [folder filename '.mat'];%[folder filename];% filename of raw datas
    disp(fullname);
    a = 1; % number of boxes to crop on the side
    w = W; % size of the last window used during PIV process
    %Dt = 1; % step between two frames that were compared during the PIV algorithm 
    %N = 0; % total number of frames processed
    %i0 = 0; % first index of the frame processed
    %b = 1; % step between frame A and A' at which velocity is computed
    
    %% Scaling 
    facq_t = freq_acq; % Frame rate in Hz
    ft = 1/facq_t ; % factor scaling for time in sec / frame
    
    % ##########################################
    %L_x = 3840; % size of the image in pixel, along larger axis (x-axis)
    %h_drone = 296.2*0.3048; % height of the drone in meter
    %theta_x = 34.15; % semi AFOV of the drone, along x-axis, in °
    %alpha_0 = 59.8; % camera pitch angle in °
    
    %dcm = 7;
    %dpx = 1192;
    
    %facq_pix = dpx/(dcm*1e-2); % scale in pixels / meter
    %facq_x = facq_pix*2/w; % scale in box / meter
    %fx = 1/facq_x; % factor scaling in meter / box
    
    % à enlever si besoin de convertir en unités de distance
    fx = 1;
    facq_pix = 1;
    % ##########################################
    
    scale_V = (facq_t/Dt) / facq_pix; % scale of the velocity in m/s
    
    % store scales in structure m
     m.scale_V = scale_V;
     m.ft = ft;
     m.fx = fx;
    
    % directory where we can save the post-processed datas
    % save_post_directory = [base 'post_processed/Post_processed_11_02_2023/'];
    % if ~exist(save_post_directory)
    %     mkdir(save_post_directory)
    % end
    
    %% load raw_datas
    
    % raw_data_path = [directory filename];
    disp(fullname);
    disp('Loading Raw Data..');
    load(fullname,'u','v','s','p');
    disp('Raw Data Loaded');
    
    disp(['Window size = ' num2str(w) ''])
    % post-processing
    [u_filt,v_filt] = PIV_banquise_postprocessing(u,v,w,N);
    
    %post_pro_fullname = [ base ]
    
    % create a matlab structure from the post-processed data file 
    m = genere_structure_banquise(u,v,u_filt,v_filt,fx,ft,scale_V,a,p,s,w);
    
    % Store parameters used for PIV processing
    m.Dt = Dt; 
    m.i0 = i0;
    m.b = b;
    %m.h_drone = h_drone;
    %m.theta_x = theta_x;
    m.facq_pix = facq_pix;
    %m.alpha_0 = alpha_0;
    
    % save the structure under a new file name
    disp('saving the new data file, with name:');
    disp(replace(fullname,'.mat','_total_processed.mat'));
    %save_name = fullname;
    %save_name = replace(save_name,'.mat','_total_processed.mat');
    save(replace(fullname,'.mat','_total_processed.mat'),'m','-v7.3');
    disp('2nd step done.');
    
    %% Plot the demodulated field at the forcing frequency
    %%
    %fig_folder = 'X:/Banquise/Vasco/Frigo_pmmh/29032024/100Hz/video_demod/';
    demod_folder = [folder 'video_demod_W' num2str(w) '_Dt' num2str(Dt) '/'];
    if ~exist(demod_folder)
        mkdir(demod_folder)
    end
    fig_folder = [demod_folder 'img_seq/'];
    if ~exist(fig_folder)
        mkdir(fig_folder)
    end
    disp('video folder exists')
    %% Création de la projection radiale
    % definition des coordonnées de la source:
    %xs = 0 ;
    %ys = 0 ;
    %xs=0;
    %ys=0;
    xs = 60 ;
    ys = -10 ;
    [Y,X] = meshgrid(m.y,m.x);
    Ux = (1./sqrt((X-xs).^2 + (Y-ys).^2)) .* (X - xs);
    Uy = (1./sqrt((X-xs).^2 + (Y-ys).^2)) .* (Y - ys);
    Ux3d = repmat(Ux,[1 1 size(m.Vx,3)]);
    Uy3d = repmat(Uy,[1 1 size(m.Vy,3)]);
    % Le signal d'intéret est le suivant: (projection radiales des vitesses)
    %H = Ux3d .* m.Vx + Uy3d .* m.Vy;
    H = m.Vy;
    disp('H array created')
    
    %%
    demodulate = 1;
    if demodulate==1
    %load('PIV_processed_i00_N0_Dt1_b1_W16_full_total_processed');
    
    % Etape 1 : initialisations
    
    Y=zeros(size(H,1),size(H,2));
    Nframes=size(H,3);
    
    % Vecteur temps
    
    
    
    t=[0:Nframes-1]/freq_acq;
    
    phase=exp(2*sqrt(-1)*pi*f_exc*t);
    
    % Construction de la matrice de demodulation
    clear A;
    A=repmat(phase,size(H,2),1);
    B(1,:,:)=A;
    demod=repmat(B,size(H,1),1,1);
    clear A;
    clear B;
    
    % On fait ensuite le produit entre H et demod (terme à terme), et la somme
    % dans la direction temporelle pour sortir le coef de fourier
    disp('demodulation en cours');
    output=sum(H.*demod,3);
    %output2 = output-mean(output,'all');
    disp('demodulation teminée')
    end
    %% plots
    if demodulate==1
    %demod_fig = figure;
    %imagesc(real(output)./abs(output))
    %colorbar()
    %caxis([0 15000])
    
    %figname = [fig_folder 'demodulated_field'];
    
    %saveas(demod_fig,figname,'fig')
    
    N_frames=48;
    t=[0:N_frames-1]/freq_acq;
    
    nb_phases = 1;
    
    figdataname = sprintf("figdata_complex.mat");
    data = output; %data = output2;
    save(demod_folder+figdataname,'data','-v7.3');
    
        figure;
    
    for i=1:N_frames*nb_phases
        phase=exp(-2*sqrt(-1)*pi*i/N_frames);
        imagesc(m.x,flip(m.y),transpose(real(output.*phase)./abs(output)));%imagesc(m.x,flip(m.y),transpose(real(output2.*phase)./abs(output2)));
        caxis([-1 1])
        colorbar()
        daspect([1 1 1])
        figname = sprintf("figure_%d.tiff",i);%[fig_folder 'demodulated_field_' string(i)];
        saveas(gcf,fig_folder + figname)
        pause(0.01)
        i
        %T(i) = getframe(gcf);
    
    
    %     close(gcf);
    end 
    end
end
