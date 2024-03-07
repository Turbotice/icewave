
%% Definition of the folder to process and where to save the results

date = '20240223';
base_img = ['/media/turbots/Hublot24/Share_hublot/PIV_images/' date(5:end) '/Drones/mesange/'];
base = ['/media/turbots/Hublot24/Share_hublot/Data/' date(5:end) '/Drones/mesange/'];
% base = '/media/turbots/DATA/thiou/labshared2/Banquise/Rimouski_2023/Data/drone/20230310/contexte/video/';
% base = 'W:/Banquise/Rimouski_2023/Data/drone/20230310/contexte/video/';
% folder = [base date '/contexte/video/'];
folder_img = [base_img 'frac+waves_002/']; 
filelist = dir([folder_img 'im*.tiff']); dirnames={};
[dirnames{1:length(filelist),1}] = deal(filelist.name);
dirnames = sortrows(dirnames);
amount = length(dirnames);
disp(amount)
%save_folder = '/media/turbots/DATA/thiou/labshared2/Banquise/Rimouski_2023/Traitement_donnees/PIV_Sebastien';
save_folder = [base 'matData/frac+waves_002/'];
%dirsave=[save_folder 'raw_datas/'];
dirsave = save_folder;
for i=1:length(filelist)
    disp([num2str(i) ' : ' filelist(i).name])
end

%% Main script

if exist(dirsave,'dir') ~= 7
    mkdir(dirsave)
end
%

% pause(10000)
for i=1:1%amount
    folderfile = filelist(i).folder;
    name=filelist(i).name;%strcat(folder, '/', dirnames{i,1},'/cam_2/');
    disp(name)
    %    jetstring=extractBetween(dirnames{i}, 'jet_', '_f');
    %    jet=strcmp(extractBetween(dirnames{i}, 'jet_', '_f'), 'on');
    %    freq=extractBetween(name, strcat(jetstring, '_f'), 'Hz')
    %    f=str2double(strrep(freq '_', '.'))
    
    
    forcing='';%extractBetween(dirnames{i}, 'mikro50mm_', '_z');
%     namesave= extractBetween([name '.'],'_offset','.');
    namesave = dirnames{i,1};
    %namesave = num2str(i);
    %disp(forcing)
    
    %    matname=string(strcat(['D:\Surface waves\Wave jet interaction\' date '\matData\'],'jet_', jetstring,'_f_', freq,'Hz_down.mat'));
        %namesave=extractBetween(name, 'z740mm_', '_1');
        
    i0 = 150; %process starting from image i0
    N = 0; %number of frames to analyze
    Dt = 4; %ratio between the fps and the scanning frequency (number of image between image A and image B)
    b = 1; %number of images between image A and image A' (from one step to an other)
    xROI = 1 ;
    widthroi = 3840;
    yROI = 1;
    heightroi = 2160;
    W = 32;
    %namesave = [namesave 'N_' num2str(N) '_i0_' num2str(i0)];
    prefix = ['PIV_processed_i0' num2str(i0) '_Dt' num2str(Dt) '_b' num2str(b) '_W' num2str(W) '_full'];
    matname = string(strcat(dirsave,prefix,'.mat'));
    disp(matname)

    if ~isfile(matname)
        disp(name)
        fullname = [folderfile];% '/' name];
%        [u,v,s,p] = PIVlab_commandline_movie_banquise(fullname,'',N,i0,Dt,b);
        [u,v,s,p]=PIVlab_commandline_parallel_stereo(fullname,'',N,i0,Dt,b);

        %[u,v,s,p]=PIVlab_commandline(name);
        
        % save variables u , v ,s and p in a .mat file
        save(matname,'u','v','s','p')
        
    else
        disp(matname)
        disp('Mat file already exist, change matname to run')
    end
    disp(i)
end