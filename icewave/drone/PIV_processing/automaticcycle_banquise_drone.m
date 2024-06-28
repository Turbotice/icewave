
%% Definition of the folder to process and where to save the results

date = '20240223';
% base where images are saved and where we want to save data 
base_img = ['/media/turbots/DATA/thiou/labshared2/SagWin2024/Data/' date(5:end) '/Drones/mesange/'];
base = ['/media/turbots/DATA/thiou/labshared2/SagWin2024/Data/' date(5:end) '/Drones/mesange/'];

folder_img = [base_img '35-waves_014/35-waves_014/']; % folder where images are saved 
filelist = dir([folder_img 'im*.tiff']); dirnames={};
% sort images 
[dirnames{1:length(filelist),1}] = deal(filelist.name);
dirnames = sortrows(dirnames);
amount = length(dirnames);
disp(amount)
%save_folder = '/media/turbots/DATA/thiou/labshared2/Banquise/Rimouski_2023/Traitement_donnees/PIV_Sebastien';
dirsave = [base 'matData/35-waves_014/'];
for i=1:length(filelist)
    disp([num2str(i) ' : ' filelist(i).name])
end

%% Main script

if exist(dirsave,'dir') ~= 7
    mkdir(dirsave)
end

for i=1:1%amount
    folderfile = filelist(i).folder;
    name=filelist(i).name;
    disp(name)
        
    forcing='';%extractBetween(dirnames{i}, 'mikro50mm_', '_z');
%     namesave= extractBetween([name '.'],'_offset','.');
    namesave = dirnames{i,1};
    %namesave = num2str(i);
    %disp(forcing)
    
    % Define parameters to process PIV
    i0 = 0; %process starting from image i0
    N = 0; %number of frames to analyze
    Dt = 4; %ratio between the fps and the scanning frequency (number of image between image A and image B)
    b = 1; %number of images between image A and image A' (from one step to an other)
    ROI.x = 650 ;
    ROI.width = 3190;
    ROI.y = 1;
    ROI.height = 2159;
    W = 32;
    prefix = ['PIV_processed_i0' num2str(i0) '_Dt' num2str(Dt) '_b' num2str(b) '_W' num2str(W) ...
        '_xROI' num2str(ROI.x) '_width' num2str(ROI.width) '_yROI' ...
        num2str(ROI.y) '_height' num2str(ROI.height)];
    matname = string(strcat(dirsave,prefix,'.mat'));
    disp(matname)

    if ~isfile(matname)
        disp(name)
%        [u,v,s,p] = PIVlab_commandline_movie_banquise(fullname,'',N,i0,Dt,b);
        [u,v,s,p]=PIVlab_commandline_parallel_stereo(folderfile,'',N,i0,Dt,b,ROI);

        % save variables u , v ,s and p in a .mat file
        save(matname,'u','v','s','p')
        
    else
        disp(matname)
        disp('Mat file already exist, change matname to run')
    end
    disp(i)
end
