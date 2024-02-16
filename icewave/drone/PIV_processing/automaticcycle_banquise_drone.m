
%% Definition of the folder to process and where to save the results

date = '20230310';
base = '/media/turbots/DATA/thiou/labshared1/Banquise/Sebastien/DJI_0308/';
% base = '/media/turbots/DATA/thiou/labshared2/Banquise/Rimouski_2023/Data/drone/20230310/contexte/video/';
% base = 'W:/Banquise/Rimouski_2023/Data/drone/20230310/contexte/video/';
% folder = [base date '/contexte/video/'];
folder = base;
filelist = dir([folder '*_images']); dirnames={};
[dirnames{1:length(filelist),1}] = deal(filelist.name);
dirnames = sortrows(dirnames);
amount = length(dirnames);
disp(amount)
%save_folder = '/media/turbots/DATA/thiou/labshared2/Banquise/Rimouski_2023/Traitement_donnees/PIV_Sebastien';
save_folder = '/media/turbots/DATA/thiou/labshared1/Banquise/Sebastien/Traitement_drone_20230310/matData/raw_datas_DJI_0308/';
% save_folder = 'Y:/Banquise/Sebastien/Traitement_drone_20230310/matData/raw_datas_DJI_0402/';
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
        
    i0 = 0; %process starting from image i0
    N = 0; %process the entire number of images in the tiff folder
    Dt = 4; %ratio between the fps and the scanning frequency (number of image between image A and image B)
    b = 1; %number of images between image A and image A'
    xROI = 1 ;
    widthroi = 3840;
    yROI = 1;
    heightroi = 2160;
    W = 64;
    %namesave = [namesave 'N_' num2str(N) '_i0_' num2str(i0)];
    prefix = ['PIV_processed_Dt' num2str(Dt) '_b' num2str(b) '_W' num2str(W) '_full_'];
    matname=string(strcat(dirsave,prefix, namesave,'.mat'));
    disp(matname)
    if ~isfile(matname)
        disp(name)
        fullname = [folderfile '/' name];
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