
%% Definition of the folder to process and where to save the results

date = '0211';
drone_name = 'mesange';
exp_ID = '08-stereo_002';
% base where images are saved and where we want to save data 

base_img = ['/media/turbots/Backup25/Data/PIV_images/' date '/Drones/' drone_name '/'];
base_save = ['/media/turbots/Backup25/Data/' date '/Drones/' drone_name '/'];

folder_img = [base_img  exp_ID '/']; % folder where images are saved 
filelist = dir([folder_img 'im*.tiff']); dirnames={};
% sort images 
[dirnames{1:length(filelist),1}] = deal(filelist.name);
dirnames = sortrows(dirnames);
amount = length(dirnames);
disp(amount)

dirsave = [base_save 'matData/' exp_ID '/'];
for i=1:length(filelist)
    disp([num2str(i) ' : ' filelist(i).name])
end

%% Main script

if exist(dirsave,'dir') ~= 7
    mkdir(dirsave)
end

% Define parameters to process PIV
i0 = 0; %process starting from image i0
N = 0; %last frame to analyze
Dt = 7; %ratio between the fps and the scanning frequency (number of image between image A and image B)
b = 1; %number of images between image A and image A' (from one step to an other)
ROI.x = 1 ;
ROI.width = 3839;
ROI.y = 1;
ROI.height = 2159;
w = 32; 

if w == 32
    nstep = 3;
elseif w == 64
    nstep = 2;
end 

% Standard PIV Settings
s = cell(15,2); % To make it more readable, let's create a "settings table"
%Parameter                          %Setting           %Options
s{1,1}= 'Int. area 1';              s{1,2}=128;         % window size of first pass
s{2,1}= 'Step size 1';              s{2,2}=32;         % step of first pass
s{3,1}= 'Subpix. finder';           s{3,2}=2;          % 1 = 3point Gauss, 2 = 2D Gauss
s{4,1}= 'Mask';                     s{4,2}=[];         % If needed, generate via: imagesc(image); [temp,Mask{1,1},Mask{1,2}]=roipoly;
s{5,1}= 'ROI';                      s{5,2}=[ROI.x,ROI.y,ROI.width,ROI.height];         % Region of interest: [x,y,width,height] in pixels, may be left empty
%s{5,1}= 'ROI';                      s{5,2}=[];         % Region of interest: [x,y,width,height] in pixels, may be left empty
s{6,1}= 'Nr. of passes';            s{6,2}= nstep;          % 1-4 nr. of passes. Each path is achieved with a specific interrogation area
s{7,1}= 'Int. area 2';              s{7,2}=64;        % second pass window size
s{8,1}= 'Int. area 3';              s{8,2}=32;         % third pass window size
s{9,1}= 'Int. area 4';              s{9,2}=32;         % fourth pass window size
s{10,1}='Window deformation';       s{10,2}='*spline'; % '*spline' is more accurate, but slower
s{11,1}='Repeated Correlation';     s{11,2}=0;         % 0 or 1 : Repeat the correlation four times and multiply the correlation matrices.
s{12,1}='Disable Autocorrelation';  s{12,2}=0;         % 0 or 1 : Disable Autocorrelation in the first pass.
s{13,1}='Correlation style';        s{13,2}=0;         % 0 or 1 : Use circular correlation (0) or linear correlation (1).
s{14,1}='Repeat last pass';         s{14,2}=0; % 0 or 1 : Repeat the last pass of a multipass analyis
s{15,1}='Last pass quality slope';   s{15,2}=0.025; % Repetitions of last pass will stop when the average difference to the previous pass is less than this number.

% Standard image preprocessing settings
p = cell(10,1);
%Parameter                       %Setting           %Options
p{1,1}= 'ROI';                   p{1,2}=s{5,2};     % same as in PIV settings
p{2,1}= 'CLAHE';                 p{2,2}=1;          % 1 = enable CLAHE (contrast enhancement), 0 = disable
p{3,1}= 'CLAHE size';            p{3,2}=64;         % CLAHE window size
p{4,1}= 'Highpass';              p{4,2}=0;          % 1 = enable highpass, 0 = disable
p{5,1}= 'Highpass size';         p{5,2}=15;         % highpass size
p{6,1}= 'Clipping';              p{6,2}=0;          % 1 = enable clipping, 0 = disable
p{7,1}= 'Wiener';                p{7,2}=0;          % 1 = enable Wiener2 adaptive denoise filter, 0 = disable
p{8,1}= 'Wiener size';           p{8,2}=3;          % Wiener2 window size
p{9,1}= 'Minimum intensity';     p{9,2}=0.0;        % Minimum intensity of input image (0 = no change)
p{10,1}='Maximum intensity';     p{10,2}=1.0;       % Maximum intensity on input image (1 = no change)

% create a structure for raw results 
PIV_param = struct('i0',i0,'N',N,'Dt',Dt,'b',b,'ROI',ROI,'w',w);
PIV_param.p = p;
PIV_param.s = s;

for i=1:1%amount
    folderfile = filelist(i).folder;
    name=filelist(i).name;
    disp(name)
        
    forcing='';%extractBetween(dirnames{i}, 'mikro50mm_', '_z');

    namesave = dirnames{i,1};
    
    prefix = ['PIV_processed_i0' num2str(i0) '_N' num2str(N) '_Dt' num2str(Dt) '_b' num2str(b) '_W' num2str(w) ...
        '_xROI' num2str(ROI.x) '_width' num2str(ROI.width) '_yROI' ...
        num2str(ROI.y) '_height' num2str(ROI.height)];
    matname = string(strcat(dirsave,prefix,'.mat'));
    disp(matname)

    if ~isfile(matname)
        disp(name)
        
        [u,v,x_pix,y_pix]=PIVlab_commandline_parallel_stereo(folderfile,'',N,i0,Dt,b,s,p);

            
        % save variables u , v , x and y pixel coordinate of boxes as well as
        % PIV parameters in a .mat file
        save(matname,'u','v','x_pix','y_pix','PIV_param')
        
    else
        disp(matname)
        disp('Mat file already exist, change matname to run')
    end
    disp(i)
end
