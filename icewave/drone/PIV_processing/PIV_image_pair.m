% Perform PIV on a pair of images 

% directory where images are saved 
directory = 'W:/SagWin2024/Data/0223/Drones/bernache/12-waves_010/12-waves_010/'; 

% name of the two images 
filename1 = 'im_0000.tiff';
filename2 = 'im_0002.tiff';

% define ROI structure 
ROI.x = 1 ;
ROI.width = 3839;
ROI.y = 1;
ROI.height = 2159;

% Standard PIV Settings
s = cell(11,2); % To make it more readable, let's create a "settings table"
%Parameter                          %Setting           %Options
s{1,1}= 'Int. area 1';              s{1,2}=128;         % window size of first pass
s{2,1}= 'Step size 1';              s{2,2}=64;         % step of first pass
s{3,1}= 'Subpix. finder';           s{3,2}=2;          % 1 = 3point Gauss, 2 = 2D Gauss
s{4,1}= 'Mask';                     s{4,2}=[];         % If needed, generate via: imagesc(image); [temp,Mask{1,1},Mask{1,2}]=roipoly;
s{5,1}= 'ROI';                      s{5,2}=[ROI.x,ROI.y,ROI.width,ROI.height];         % Region of interest: [x,y,width,height] in pixels, may be left empty
%s{5,1}= 'ROI';                      s{5,2}=[];         % Region of interest: [x,y,width,height] in pixels, may be left empty
s{6,1}= 'Nr. of passes';            s{6,2}=3;          % 1-4 nr. of passes. Each path is achieved with a specific interrogation area
s{7,1}= 'Int. area 2';              s{7,2}=64;        % second pass window size
s{8,1}= 'Int. area 3';              s{8,2}=32;         % third pass window size
s{9,1}= 'Int. area 4';              s{9,2}=32;         % fourth pass window size
s{10,1}='Window deformation';       s{10,2}='*spline'; % '*spline' is more accurate, but slower
s{11,1}='Repeated Correlation';     s{11,2}=0;         % 0 or 1 : Repeat the correlation four times and multiply the correlation matrices.
s{12,1}='Disable Autocorrelation';  s{12,2}=0;         % 0 or 1 : Disable Autocorrelation in the first pass.
s{13,1}='Correlation style';        s{13,2}=0;         % 0 or 1 : Use circular correlation (0) or linear correlation (1).
s{14,1}='Repeat last pass';   s{14,2}=0; % 0 or 1 : Repeat the last pass of a multipass analyis
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

% Number of cores 
nr_of_cores = 1;

% Perform PIV analysis 
[x0, y0, u0, v0, typevector,correlation_map] = ...
			piv_analysis(directory, filename1, filename2,p,s,nr_of_cores,false);



