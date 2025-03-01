%% Script aims at georeferencing images extracted from a drone video of format 4K
package_path = 'C:/Users/sebas/icewave/icewave/drone/Function_package/';
addpath(package_path, '-end');


year = '2024';
date = '0921';
drone_ID = 'bernache';
path2data =  ['C:/Users/sebas/Desktop/Amundsen_RA_2024/Data/' year '/' date '/Drones/' drone_ID];

filelist = dir([path2data '/' '*ortho_003' ]);
img_format = '*.tiff';
feet2meter = 1/3.2808; % 1 feet -> x meter

% fields to be converted
keys = {'h','alpha'}; % keys for which unit conversion is needed
initial_units = {'ft','deg'}; % initial units of drone parameters
IS_units = {'m','rad'}; % IS units of drone parameters
coef_conversion = [feet2meter pi/180]; % conversion coefficients from initial unit to IS units

% Create a folder where results can be saved
fig_folder = [filelist(1).folder '\' filelist(1).name '\Figures\'];
if ~isfolder(fig_folder)
    mkdir(fig_folder)
end
fig_resolution = 300; % resolution of saved fig


% Parameters for image binarization
default_sensitivity = 0.58; % sensitivity used by adaptthresh 
radius_open = 3; % radius of structuring element to perform erosion
minArea = 100; % minimal area (in pixels) of detected objects
conn = 4; % connection used for bwboundaries


for folder_idx = 1 : 1
    current_folder = filelist(folder_idx).name;
    path_folder = [filelist(folder_idx).folder '\' current_folder];

    % load drone parameters in a structure 
    filelist_drone_param = dir([path_folder '\' '*.txt']);
    for file_idx = 1 : length(filelist_drone_param)
        current_filename = filelist_drone_param(file_idx).name;
        lines = readcell([filelist_drone_param(file_idx).folder '\' current_filename]);
        s = struct();
        s.units = struct();
        for variable_idx = 1: size(lines,1)
            header = lines{variable_idx,1};
            s.(header) = lines{variable_idx,3};
            s.units.(header) = lines{variable_idx,4};
        end 

        % Change drone parameters coefficient 
        s = conversion_structure_units(s,keys,initial_units,IS_units,coef_conversion);

        % load JPG images one by one 
        filelist_image = dir([path_folder '\selection_frames\' img_format]);
        for idx = 1 : length(filelist_image)
            img_name = [filelist_image(idx).folder '\' filelist_image(idx).name];
            disp(img_name)

            current_img = imread(img_name);
            
            txt_cell = split(filelist_image(idx).name,'.');
            prefix_current_fig = txt_cell{1};
           
            % BW_param = struct('sensitivity',default_sensitivity,'radius_open',radius_open,...
            %     'minArea',minArea,'conn',conn);
            % 
            % results = ice_floes_detection(current_img,s,BW_param,fig_folder,prefix_current_fig,fig_resolution);
            % 
            % results.figfilename = img_name;
            % 
            % file2load = [fig_folder prefix_current_fig '_Data_sensitivity' num2str(results.BW_param.sensitivity) ...
            %     '_imopen' num2str(results.BW_param.radius_open) '_conn' num2str(results.BW_param.conn) '_minarea' ...
            %     num2str(results.BW_param.minArea) '.mat'];
            % 
            % disp('Saving Data..')
            % save(file2load,'results','-v7.3')
            % disp('Data saved !')
           
        end 

    end 

end 



%% Load Data

file2load = [fig_folder prefix_current_fig '_Data_sensitivity' num2str(results.BW_param.sensitivity) ...
    '_imopen' num2str(results.BW_param.radius_open) '_conn' num2str(results.BW_param.conn) '_minarea' ...
    num2str(results.BW_param.minArea) '.mat'];

disp('Loading Data..')
load(file2load)
disp('Data loaded !')

%% Create histogram from data

A = [];
stats = results.stats;
for k = 1 : length(stats)
    A = cat(1,A,stats(k).area_real);
end 
Amin_real = 100; % minimal area of ice floes in real space
Amax_real = max(A);

Nbins = 30;
figure(7),
h = histogram(A,Nbins,'BinLimits',[Amin_real Amax_real]);
xlabel('Floe Area $\rm (m^{2})$')
ylabel('Counts')
title(['FSD - $\rm{Area} >' num2str(Amin_real) ' \: \rm m^2$'])
set(findall(gcf,'-property','FontSize'),'FontSize',15)
set_Papermode(gcf)

figname_cell = split(results.figfilename,'\');
prefixe = figname_cell{end}(1:6);
disp(prefixe)

figname = [fig_folder prefixe '_Histogram_Amin_real' num2str(Amin_real) '_Nbins' num2str(Nbins)];
saveas(gcf,[figname '.fig'])
exportgraphics(gcf,[figname '.png'],'Resolution',fig_resolution)
exportgraphics(gcf,[figname '.pdf'],'Resolution',fig_resolution)

% #######################################
%% MAIN DEVELOPMENTS 
% #######################################

% ############################################
%% Cut picture into smaller subpicture
% ############################################

figure(1)
imshow(current_img)
title('Original RGB image')

img_gray = rgb2gray(current_img);

[Ly,Lx,~] = size(current_img);
nx = 2; % number of sub image in x-direction  
ny = 2; % number of sub image in y-direction 

default_sensitivity = 0.6; % default sensitivity used for adapt thresholding 
default_disk = 3; % default disk radius for image opening
default_minArea = 50; % default minimal area of detected objects on picture 
default_h = 4; % default value for extended-minimal transformation 

sub_Lx = floor(Lx/nx);
sub_Ly = floor(Ly/ny);
array_dim_1 = int16(ones(ny,1) .* sub_Ly);
array_dim_2 = int16(ones(nx,1) .* sub_Lx);

if (rem(Lx,nx) > 0) || (rem(Ly,ny) > 0)

    array_dim_1(end) = array_dim_1(end) + rem(Ly,ny);
    array_dim_2(end) = array_dim_2(end) + rem(Lx,nx);

end 

% converts the image into 4 different matrices 
C = mat2cell(current_img,array_dim_1 , array_dim_2, 3);
matrix_BW = cell(size(C));



%%
% select first sub image 
for sub_idx = 1:numel(C)

    sub_img = C{sub_idx};
    
    figure(2)
    imshow(sub_img)
    title('Sub image')
    
    % select reference RGB values
    disp('Select points on the figure to set reference color')
    [x,y] = getpts(gcf); 
    x = int16(x);
    y = int16(y);
    c = zeros([size(x,1),3]); % array of colors
    for k = 1: length(x)
        disp(['RGB values : R = ' num2str(sub_img(y(k),x(k),1)) ', G = ' num2str(sub_img(y(k),x(k),2))...
            ', B = ' num2str(sub_img(y(k),x(k),3))])
    
        c(k,:) = sub_img(y(k),x(k),:);
    end 
    
    c0 = mean(c,1);
    disp(['Averaged RGB values of selected points : ' num2str(c0)])
    
    c0 = c0 / 255 ; 
    img_double = im2double(sub_img); % convert image matrix to double
    I = sqrt((img_double(:,:,1) - c0(1)).^2 + (img_double(:,:,2) - c0(2)).^2 + (img_double(:,:,3) - c0(3)).^2);
    I = ones(size(I)) - I;
    
    figure(3), 
    imagesc(I)
    colormap("gray")
    colorbar()
    title('Distance to reference color')
    
    T = graythresh(I);
    BW = imbinarize(I,T);
    figure(4),
    imshow(BW)
    title('Binarized sub-image - Otsu method')

    prompt = 'Do you want to use an adapted thresholding ? y/n [n]';
    txt_BW = input(prompt,"s");

    if strcmp(txt_BW,'y')
    
        % Binarize using adapt threshold
        test_sensitivity = 'n';
        while strcmp(test_sensitivity,'n')
            prompt = ['Set value of sensitivity (default ' num2str(default_sensitivity) ') : '];
            txt = input(prompt,"s");
            if isempty(txt)
                sensitivity = default_sensitivity;  % sensitivity used for adaptthreshold
            else 
                sensitivity = str2double(txt);
            end 
            disp(['Current sensitivity = ' num2str(sensitivity)])
        
    
            T = adaptthresh(I,sensitivity,'ForegroundPolarity','bright','NeighborhoodSize',2*floor(size(I)/16) + 1);
            BW = imbinarize(I,T);
        
            imshow(BW)
            title('Binarized sub-image')
    
            prompt = 'Is sensitivity sufficient ? y/n [y]';
            txt = input(prompt,"s");
            if isempty(txt)
                test_sensitivity = 'y';
            elseif strcmp(txt,'y')
                test_sensitivity = 'y';
            end 
        
        end 

    end

    % Image opening
    test_opening = 'n';
    while strcmp(test_opening,'n')
        prompt = ['Set value of disk radius for opening operation (default ' num2str(default_disk) ') : '];
        txt = input(prompt,"s");
        if isempty(txt)
            r_disk = default_disk;  % sensitivity used for adaptthreshold
        else 
            r_disk = str2double(txt);
        end 
        disp(['Current opening disk : ' num2str(r_disk)])
        se = strel('disk',r_disk);
        Io = imopen(BW,se);
        
        Ioverlay = labeloverlay(I,Io);
        figure(5)
        imshow(Ioverlay)
        title('Opening binarized sub-image')

        prompt = 'Is opening sufficient ? y/n [y]';
        txt = input(prompt,"s");
        if isempty(txt)
            test_opening = 'y';
        elseif strcmp(txt,'y')
            test_opening = 'y';
        end
    end 

    % filter small objects
    test_minarea = 'n';
    while strcmp(test_minarea,'n')
        prompt = ['Set value of minimal area od detected objects (default : ' num2str(default_minArea) ') : '];
        txt = input(prompt,"s");
        if isempty(txt)
            minArea = default_minArea;  % sensitivity used for adaptthreshold
        else 
            minArea = str2double(txt);
        end 

        Io_filt = bwareaopen(Io,minArea);
        Ioverlay = labeloverlay(I,Io_filt);

        figure(6)
        imshow(Ioverlay)
        title('Binarized image after filtering')

        prompt = 'Is filtering sufficient ? y/n [y]';
        txt = input(prompt,"s");
        if isempty(txt)
            test_minarea = 'y';
        elseif strcmp(txt,'y')
            test_minarea = 'y';
        end

    end

    % Store binarized image in a cell 
    matrix_BW{sub_idx} = Io_filt;

end 

%% Build full BW image 

img_gray = rgb2gray(current_img);
full_BW = img_gray;
for i = 1 : ny
    for j = 1 : nx
        
    full_BW(1 + (i-1)*sub_Ly:sub_Ly + (i-1) * sub_Ly, 1 + (j-1)*sub_Lx : sub_Lx + (j-1)*sub_Lx) ...
        = matrix_BW{i,j};

    end 
end 

full_overlay = labeloverlay(img_gray,full_BW);
figure, 
imshow(full_overlay)


%% Detect objects using bwboundaries 
conn = 4;
[B,Lbw,n,~] = bwboundaries(full_BW,conn,'noholes','CoordinateOrder',"xy");
disp(['Number of detected objects :' num2str(n)])

statsbw = regionprops(Lbw,'Centroid','BoundingBox','PixelList','PixelIdxList');

figure(7), 
imshow(img_gray)
hold on 
for k = 1 : length(statsbw)

    boundary = B{k};
    statsbw(k).('boundary_pix') = boundary;
    plot(boundary(:,1),boundary(:,2),'r')
    hold on 

end 

Lrgb = label2rgb(Lbw,"jet","w","shuffle");

figure(7);
imshow(img_gray)
hold on
himage = imshow(Lrgb);
himage.AlphaData = 0.3;
title("Colored Labels Superimposed Transparently on Original Image")



%% Zoom on a single ice floe

main_fig = figure(8);
imshow(img_gray)
hold on 
for k = 1 : length(statsbw)

    plot(statsbw(k).boundary_pix(:,1),statsbw(k).boundary_pix(:,2),'r')
    hold on 
end

sz = size(full_BW);

manual_segmentation = 1;
while manual_segmentation
    prompt = 'Do you want to perform manual segmentation ? y/n [y]';
    txt = input(prompt,"s");
    if isempty(txt)
        manual_segmentation = 1;  % sensitivity used for adaptthreshold
    elseif strcmp(txt,'n')
        manual_segmentation = 0;
    end 

    disp('Select undersegmented ice floes')
    [xfloe,yfloe] = getpts(main_fig);
    linear_idx = sub2ind(sz,yfloe,xfloe);
    corresponding_floe = find_objects(xfloe,yfloe,statsbw,sz);
    
    hold on 
    plot(statsbw(corresponding_floe).Centroid(1),statsbw(corresponding_floe).Centroid(2),'ro')

    % Get BoundingBox of selected ice floe
    
    window = [statsbw(corresponding_floe).BoundingBox(2) ...
        statsbw(corresponding_floe).BoundingBox(2) + statsbw(corresponding_floe).BoundingBox(4) ; ... 
        statsbw(corresponding_floe).BoundingBox(1) ...
        statsbw(corresponding_floe).BoundingBox(1) + statsbw(corresponding_floe).BoundingBox(3) ];
    
    window = round(window);
    
    
    sub_img = img_gray(window(1,1) : window(1,2),window(2,1) : window(2,2));
    sub_bw = full_BW(window(1,1) : window(1,2),window(2,1) : window(2,2));
    label_img = labeloverlay(sub_img,sub_bw,'Transparency',0.8);
    
    xc = statsbw(corresponding_floe).Centroid(1) - statsbw(corresponding_floe).BoundingBox(1);
    yc = statsbw(corresponding_floe).Centroid(2) - statsbw(corresponding_floe).BoundingBox(2);
    % figure(9)
    % imshow(label_img)
    % hold on 
    % scatter(xc,yc,'ro')
    % close(9)
    % Split selected ice floe 
    
    D = -bwdist(~sub_bw);
    
    % Use imextendedmin
    test_extendedmin = 'n';
    while strcmp(test_extendedmin,'n')
            prompt = ['Set value for extended-minima transform (default ' num2str(default_h) ') : '];
            txt = input(prompt,"s");
            if isempty(txt)
                h = default_h;  % sensitivity used for adaptthreshold
            else 
                h = str2double(txt);
            end 
            disp(['Current h : ' num2str(h)])
            
            mask = imextendedmin(D,h);
            figure(10)
            imshow(labeloverlay(sub_img,mask,'Transparency',0.5));
    
            prompt = 'Is extended-minima transform correct ? y/n [y]';
            txt = input(prompt,"s");
            if isempty(txt)
                test_extendedmin = 'y';
            elseif strcmp(txt,'y')
                test_extendedmin = 'y';
            end
    end

    D2 = imimposemin(D,mask);

    % Apply watershed algorithm 
    
    subL = watershed(D2);
    subL(~sub_bw) = 0;
    
    sub_stats = regionprops(subL,'Centroid','BoundingBox','PixelList','PixelIdxList');
    subLrgb = label2rgb(subL,"jet","w","shuffle");
    
    figure(11),
    imshow(sub_img)
    hold on
    himage = imshow(subLrgb);
    himage.AlphaData = 0.3;
    title("Colored Labels Superimposed Transparently on Original Image")

    disp('Select new ice floes')
    [xfloe,yfloe] = getpts(gcf);
    new_floes_idx = find_objects(xfloe,yfloe,sub_stats,size(sub_img));
    close(11)

    % get boundaries of respective new ice floes 
    mask_bound = zeros(size(subL));
    for k = 1 : length(new_floes_idx)
        idx = new_floes_idx(k);
        mask_bound(sub_stats(idx).PixelIdxList) = 1;
    end 
    
    new_boundaries = bwboundaries(mask_bound,conn,'noholes','CoordinateOrder',"xy");
    % Add boundaries of new objects
    for i0 = 1 : length(new_floes_idx)
        idx = new_floes_idx(i0);
        sub_stats(idx).boundary_pix = new_boundaries{i0};
    end 

    % replace existing floe by new floe idx 1
    % shift to apply to pixels location of new objects
    shift = [statsbw(corresponding_floe).BoundingBox(1) statsbw(corresponding_floe).BoundingBox(2)];
    first_floe = new_floes_idx(1);
    
    new_object = subimg2main_obj_propreties(sub_stats(first_floe),shift,sz);
    
    statsbw(corresponding_floe) = new_object;
    
    N = length(statsbw);
    sz = size(img_gray);
    for k = 2 : length(new_floes_idx)
        current_floe = new_floes_idx(k);
        current_object = subimg2main_obj_propreties(sub_stats(current_floe),shift,sz);
        statsbw(N + k - 1) = current_object;
    end

    disp(length(statsbw))

    % Save data
    disp('Saving data..')
    filename = [fig_folder prefix_current_fig '_Data_BW_floes_segmentation' ];
    save(filename,'statsbw','-v7.3')
    disp('Data saved')

    main_fig = figure(8);
    imshow(img_gray)
    hold on 
    for k = 1 : length(statsbw)
    
        plot(statsbw(k).boundary_pix(:,1),statsbw(k).boundary_pix(:,2),'r')
        hold on 
    end



end



% Add a break of the while loop if we do not want to perform segmentation
% with watershed 

























% ###############################################
%% Use watershed on fully reconstructed BW image 
% ###############################################

% Compute distance map
D = -bwdist(~full_BW);
figure(7),
imagesc(D)
colorbar()
clim([-20 0])

% Kills small water ridges
conn = 4;
h = 8;
D_filt = imhmin(D,h,conn);
% D_filt = D;
% D_filt(D < 0) = D_filt(D <0) + h;
% Watershed
L = watershed(D_filt);
L(~full_BW) = 0;

disp([num2str(max(L,[],'all')) ' floes detected'])

Lrgb = label2rgb(L,"jet","w","shuffle");

figure(8),
imshow(img_gray)
hold on
himage = imshow(Lrgb);
himage.AlphaData = 0.3;
title("Colored Labels Superimposed Transparently on Original Image")

%%
[Gmag,Gdir] = imgradient(L); % boundaries of detected objects 
stats = regionprops(L,'Centroid','BoundingBox','PixelList','PixelIdxList');  

sz = size(L); 

figure,
imshow(img_gray)
hold on 
for k = 1 : length(stats)

    plot(stats(k).Centroid(1),stats(k).Centroid(2),'r.')
    
    % keep only boundaries 

    BorderIdxList = stats(k).PixelIdxList(Gmag(stats(k).PixelIdxList) > 0); % linear indices of pixels belonging to the object border
    
    [row , col] = ind2sub(sz,BorderIdxList);
    
    hold on 
    scatter(col,row,10,'b.')
    % pause()
    % hold off
end 














% #########################################################################
%% Ice floes segmentation using watershed method
% step 1 : binarization using adaptthresh
% step 2 : open by reconstruction (image opening, keeping object shape)
% step 3 : compute bwdistance (energy potential)
% step 4 : watershed : L
% step 5 : regionprops
% step 6 : computes imgradient de watershed (L)
% step 7 : extract boundaries of each object

% step 8 (optional) : get rid off objects touching the image border
% #########################################################################

idx = 1;

img_name = [filelist_image(idx).folder '\' filelist_image(idx).name];
disp(img_name)

current_img = imread(img_name);

txt_cell = split(filelist_image(idx).name,'.');
prefix_current_fig = txt_cell{1};

img_gray = rgb2gray(current_img);

sensitivity = 0.58;
figure(1),
T = adaptthresh(img_gray,sensitivity,'ForegroundPolarity','bright','NeighborhoodSize',2*floor(size(img_gray)/16) + 1);
BW = imbinarize(img_gray,T);
imshowpair(img_gray,BW,'montage')

% Compare imopen and imopen / reconstruct 
r_disk = 5;

se = strel('disk',r_disk);

Ie = imerode(BW,se);
Iobr = imreconstruct(Ie,BW);
Io = imopen(BW,strel('disk',3));

figure,
imshow(Iobr)
title(['Opening-by-Reconstruction - r = ' num2str(r_disk)])

figure, 
imshow(Io)
title('Image opening')



%%
D = -bwdist(~Iobr);
figure(3),
imagesc(D)
colormap("gray")
colorbar()

%%

conn = 4;
h = 5;
for i = 1 : length(h)
    Ifilt = imhmin(D,h(i),conn);
    L = watershed(Ifilt);
    L(~Io) = 0;
    
    disp([num2str(max(L,[],'all')) ' floes detected'])

    stats = regionprops(L,'Centroid','BoundingBox','PixelList','PixelIdxList');

    disp(['Number of ice floes detected : ' num2str(length(stats))])
    figure,
    imshow(Iobr)
    hold on 
    for k = 1 : length(stats)
        plot(stats(k).Centroid(1),stats(k).Centroid(2),'r.')
    end 
    set_Papermode(gcf)
end


%%
[Gmag,Gdir] = imgradient(L); % boundaries of detected objects 
figure, 
imshow(Gmag)
    
sz = size(Iobr); 

figure,
imshow(current_img)
hold on 
for k = 1 : length(stats)

    plot(stats(k).Centroid(1),stats(k).Centroid(2),'r.')
    
    % keep only boundaries 

    BorderIdxList = stats(k).PixelIdxList(Gmag(stats(k).PixelIdxList) > 0); % linear indices of pixels belonging to the object border
    
    [row , col] = ind2sub(sz,BorderIdxList);
    
    hold on 
    scatter(col,row,20,'b.')
    % pause()
    % hold off
end 

% #########################################
%% Build my own process using watershed
% #########################################

I = rgb2gray(current_img); % gray scaled image 
figure(1),
imshow(I)
title('Gray scaled image')


% open-by-reconstruction of the grayscaled image 
r_disk = 5;
se = strel('disk',r_disk);
figure(2),
Ie = imerode(I,se);
Iobr = imreconstruct(Ie,I);
imshow(Iobr)
title("Opening-by-Reconstruction")

% Binarize using adapttrshold
sensitivity = 0.58;
T = adaptthresh(I,sensitivity,'ForegroundPolarity','bright','NeighborhoodSize',2*floor(size(img_gray)/16) + 1);
BW = imbinarize(I,T);
figure(3),
imshow(BW)
title('Thresholded image (adaptthresh)')

% Opening of the binarized image
r_disk = 5;
se = strel('disk',r_disk);
BWo = imopen(BW,se);
figure(4)
imshow(BWo)
title('Opening of BW')

% Get rid off small objects
BW_filt = bwareaopen(BWo,50);
figure(5),
imshow(BW_filt)
title('Binarized image filtered')

% Compute distance map
D = -bwdist(~BW_filt);
figure(6),
imagesc(D)
colorbar()

% Kills small water ridges
conn = 4;
h = 10;
D_filt = imhmin(D,h,conn);


% Watershed
L = watershed(D_filt);
L(~BW_filt) = 0;

disp([num2str(max(L,[],'all')) ' floes detected'])

Lrgb = label2rgb(L,"jet","w","shuffle");

figure(7),
imshow(I)
hold on
himage = imshow(Lrgb);
himage.AlphaData = 0.3;
title("Colored Labels Superimposed Transparently on Original Image")

[Gmag,Gdir] = imgradient(L); % boundaries of detected objects 
stats = regionprops(L,'Centroid','BoundingBox','PixelList','PixelIdxList');  
sz = size(L); 

figure,
imshow(current_img)
hold on 
for k = 1 : length(stats)

    plot(stats(k).Centroid(1),stats(k).Centroid(2),'r.')
    
    % keep only boundaries 

    BorderIdxList = stats(k).PixelIdxList(Gmag(stats(k).PixelIdxList) > 0); % linear indices of pixels belonging to the object border
    
    [row , col] = ind2sub(sz,BorderIdxList);
    
    hold on 
    scatter(col,row,10,'b.')
    % pause()
    % hold off
end 


% ###############################
%% Try k-means clustering 
% ###############################
img_gray = rgb2gray(current_img);
figure, 
imshow(img_gray)


[L,Centers] = imsegkmeans(I,3);
B = labeloverlay(I,L);
figure,
imshow(B)
title("Labeled Image")



% ###############################################
%% Marker Controlled watershed segmentation 
% ###############################################

idx = 1;

img_name = [filelist_image(idx).folder '\' filelist_image(idx).name];
disp(img_name)

current_img = imread(img_name);

txt_cell = split(filelist_image(idx).name,'.');
prefix_current_fig = txt_cell{1};

I = rgb2gray(current_img);

% Compute gradient magnitude
gmag = imgradient(I);
figure(1),
imshow(gmag,[])
title("Gradient Magnitude")

% Get Foreground markers
se = strel("disk",10);
Io = imopen(I,se);
% figure(2),
% imshow(Io)
% title("Opening")

figure(3),
Ie = imerode(I,se);
Iobr = imreconstruct(Ie,I);
imshow(Iobr)
title("Opening-by-Reconstruction")

% figure(4),
% Ioc = imclose(Io,se);
% imshow(Ioc)
% title("Opening-Closing")

Iobrd = imdilate(Iobr,se);
Iobrcbr = imreconstruct(imcomplement(Iobrd),imcomplement(Iobr));
Iobrcbr = imcomplement(Iobrcbr);
figure(5),
imshow(Iobrcbr)
title("Opening-Closing by Reconstruction")

fgm = imregionalmax(Iobrcbr);
I2 = labeloverlay(I,fgm);
figure(6),
imshow(I2)
title("Regional Maxima Superimposed on Original Image")

se2 = strel(ones(10,10));
fgm2 = imerode(fgm,se2);
fgm3 = imclose(fgm2,se2);
fgm4 = bwareaopen(fgm3,100);
I3 = labeloverlay(I,fgm4);
figure(7)
imshow(I3)
title("Modified Regional Maxima Superimposed on Original Image")

sensitivity = 0.58;
T = adaptthresh(Iobrcbr,sensitivity,'ForegroundPolarity','bright','NeighborhoodSize',2*floor(size(img_gray)/16) + 1);
bw = imbinarize(Iobrcbr,T);
figure(8)
imshow(bw)
title("Thresholded Opening-Closing by Reconstruction")


D = bwdist(bw);
DL = watershed(D);
bgm = DL == 0;
figure(9)
imshow(bgm)
title("Watershed Ridge Lines")

gmag2 = imimposemin(gmag, bgm | fgm4);
figure(11)
imagesc(bgm | fgm4)
colorbar()
L = watershed(gmag2);

Lrgb = label2rgb(L,"jet","w","shuffle");

figure(10),
imshow(I)
hold on
himage = imshow(Lrgb);
himage.AlphaData = 0.3;
title("Colored Labels Superimposed Transparently on Original Image")













% ################################################
%% imextendedmin and imimposemin using bwdist
% ################################################

full_overlay = labeloverlay(img_gray,full_BW);
figure(1), 
imshow(full_overlay)

D = -bwdist(~full_BW);
figure(2)
imshow(D,[])

mask = imextendedmin(D,20);
figure(3)
imshow(labeloverlay(img_gray,mask))

D2 = imimposemin(D,mask);
figure(4)
surf(D2')
shading interp
colorbar()
%%
Ld2 = watershed(D2);
Ld2(~full_BW) = 0;

Lrgb = label2rgb(Ld2,"jet","w","shuffle");

figure(5),
imshow(img_gray)
hold on
himage = imshow(Lrgb);
himage.AlphaData = 0.3;
title("Colored Labels Superimposed Transparently on Original Image")

















% ##################################
%% Play with watershed algorithm 
% ##################################

% get foreground markers of a sub-image 
sub_idx = 1; 

sub_img = C{sub_idx};

figure(2)
imshow(sub_img)
title('Sub image')

% select reference RGB values
disp('Select points on the figure to set reference color')
[x,y] = getpts(gcf); 
x = int16(x);
y = int16(y);
c = zeros([size(x,1),3]); % array of colors
for k = 1: length(x)
    disp(['RGB values : R = ' num2str(sub_img(y(k),x(k),1)) ', G = ' num2str(sub_img(y(k),x(k),2))...
        ', B = ' num2str(sub_img(y(k),x(k),3))])

    c(k,:) = sub_img(y(k),x(k),:);
end 

c0 = mean(c,1);
disp(['Averaged RGB values of selected points : ' num2str(c0)])

c0 = c0 / 255 ; 
img_double = im2double(sub_img); % convert image matrix to double
I = sqrt((img_double(:,:,1) - c0(1)).^2 + (img_double(:,:,2) - c0(2)).^2 + (img_double(:,:,3) - c0(3)).^2);
I = ones(size(I)) - I;

figure(3), 
imagesc(I)
colormap("gray")
colorbar()
title('Distance to reference color')

%% Get Foreground markers using BW
sub_gray = rgb2gray(sub_img);

% adapt contrast
Ieq = adapthisteq(sub_gray); % CLAHE contrast enhancement 

sensitivity = 0.6;
T = adaptthresh(Ieq,sensitivity,'ForegroundPolarity','bright','NeighborhoodSize',2*floor(size(I)/16) + 1);
BW = imbinarize(Ieq,T);

figure,
imshow(BW)

% T = adaptthresh(I,sensitivity,'ForegroundPolarity','bright','NeighborhoodSize',2*floor(size(I)/16) + 1);
% BW = imbinarize(I,T);
% 
% figure,
% imshow(BW)
% title('Binarized sub-image Color-distance')
title('Binarized sub-image CLAHE')

bw2 = imfill(BW,'holes');
figure(2),
imshow(bw2)

se = strel('disk',5);
bw3 = imopen(bw2, se);
figure(3), 
imshow(bw3)

bw4 = bwareaopen(bw3, 40);
figure(4)
imshow(bw4)

bw4_perim = bwperim(bw4);
overlay1 = imoverlay(Ieq, bw4_perim);
figure(5)
imshow(overlay1)

mask_em = imextendedmax(Ieq,15);
Ilbl = labeloverlay(Ieq,mask_em);
figure(6)
imshow(Ilbl)

mask_em = imfill(mask_em,'holes');
figure(7)
imshow(labeloverlay(Ieq,mask_em))

mask_em = imerode(mask_em, ones(5,5));
figure(8)
imshow(labeloverlay(Ieq,mask_em))

mask_em = imclose(mask_em, ones(3,3));
figure(9)
imshow(labeloverlay(Ieq,mask_em))

mask_em = bwareaopen(mask_em, 40);
indices = mask_em & bw4_perim;
mask_em(indices) = 0;

figure(10),
imshow(mask_em)

%%
overlay2 = labeloverlay(I_eq, bw4_perim | mask_em);
figure(11)
imshow(overlay2)

%%
Ieq_c = imcomplement(Ieq);
I_mod = imimposemin(Ieq_c, ~bw4 | mask_em);
figure(11)
imshow(I_mod)

L = watershed(I_mod);
Lrgb = label2rgb(L,"jet","w","shuffle");

figure(12),
imshow(Ieq)
hold on
himage = imshow(Lrgb);
himage.AlphaData = 0.3;
title("Colored Labels Superimposed Transparently on Original Image")
% 







%%


% define pixel coordinates 
Lx = size(current_img,2);
Ly = size(current_img,1);

X_pix=repmat((1:Lx),Ly,1);
Y_pix=repmat((1:Ly)',1,Lx);
% define camera pixel center 
x0 = (Lx + 1)/2;
y0 = (Ly + 1)/2;

% Compute coordinates in real space from pixel coordinates 

[Xreal,Yreal] = projection_real_space(X_pix,Y_pix,x0,y0,s.h,s.alpha,s.focale);

figure(1);
imshow(current_img)

img_gray = rgb2gray(current_img);

figure(2),
pcolor(Xreal,Yreal,img_gray)
colormap("gray")
shading interp
axis image
xlabel('$x \: \rm (m)$')
ylabel('$y \: \rm (m)$')
set(findall(gcf,'-property','FontSize'),'FontSize',15)
set_Papermode(gcf)

txt_cell = split(filelist_image(idx).name,'.');
prefix_current_fig = txt_cell{1};
figname = [fig_folder prefix_current_fig '_Oblique_correction'];
saveas(gcf,[figname '.fig'])
exportgraphics(gcf,[figname '.png'],'Resolution',fig_resolution)
exportgraphics(gcf,[figname '.pdf'],'Resolution',fig_resolution)

%% Apply adaptative threshold to grayscale image

default_sensitivity = 0.58; % sensitivity used by adaptthresh 
radius_open = 3; % radius of structuring element to perform erosion
minArea = 100; % minimal area (in pixels) of detected objects
conn = 4; % connection used for bwboundaries

test_sensitivity = 'n';
while strcmp(test_sensitivity,'n')
    prompt = ['Set value of sensitivity (default ' num2str(default_sensitivity) ') : '];
    txt = input(prompt,"s");
    if isempty(txt)
        sensitivity = default_sensitivity;  % sensitivity used for adaptthreshold
    else 
        sensitivity = str2double(txt);
    end 
    disp(['Current sensitivity = ' num2str(sensitivity)])

    figure(3),
    T = adaptthresh(img_gray,sensitivity,'ForegroundPolarity','bright','NeighborhoodSize',2*floor(size(img_gray)/16) + 1);
    BW = imbinarize(img_gray,T);
    imshowpair(img_gray,BW,'montage')

    prompt = 'Is sensitivity sufficient ? y/n [y]';
    txt = input(prompt,"s");
    if isempty(txt)
        test_sensitivity = 'y';
    elseif strcmp(txt,'y')
        test_sensitivity = 'y';
    end 

end 

% Try imopen, imclose and areaopen 

% radius_close = 3;
struct_open = strel('disk',radius_open);
% struct_close = strel('disk',radius_close);
BW_cleaned = imopen(BW,struct_open);
% Bw_cleaned = imclose(BW_cleaned,struct_close);

% Filter out small objects based on area

BW_filtered = bwareaopen(BW_cleaned, minArea);
figure(4),
% imshowpair(img_gray,BW_filtered,'montage')
imshow(BW_filtered)
set_Papermode(gcf)
figname = [fig_folder prefix_current_fig '_BW_filtered_sensitivity' num2str(sensitivity) ... 
    '_imopen' num2str(radius_open) '_conn' num2str(conn) '_minarea' num2str(minArea)];
saveas(gcf,[figname '.fig'])
exportgraphics(gcf,[figname '.png'],'Resolution',fig_resolution)
exportgraphics(gcf,[figname '.pdf'],'Resolution',fig_resolution)

%% Detect ROI and extract there position

[B,L,n,~] = bwboundaries(BW_filtered,conn,'noholes','CoordinateOrder',"xy");
disp(['Number of detected objects :' num2str(n)])

stats = regionprops(L,'Centroid','PixelList');

figure(5),
imshow(img_gray)
hold on 
for k = 1 : n
    plot(stats(k).Centroid(1),stats(k).Centroid(2),'r.')
end 
set_Papermode(gcf)

figname = [fig_folder prefix_current_fig '_Detected_floes_centroid_sensitivity' num2str(sensitivity) ... 
    '_imopen' num2str(radius_open) '_conn' num2str(conn) '_minarea' num2str(minArea)];
saveas(gcf,[figname '.fig'])
exportgraphics(gcf,[figname '.png'],'Resolution',fig_resolution)
exportgraphics(gcf,[figname '.pdf'],'Resolution',fig_resolution)

%% Georeference picture and pixels of a given object + object area computation 
figure(6),
pcolor(Xreal,Yreal,img_gray)
colormap("gray")
shading interp

for idx_object = 1 : n
    boundary = B{idx_object}; % #1 - x  and #2 - y

    % georeference object boundaries
    boundary_real = boundary;
    [boundary_real(:,1),boundary_real(:,2)] = projection_real_space(boundary(:,1),boundary(:,2),x0,y0,s.h,s.alpha,s.focale);

    hold on 
    plot(boundary_real(:,1),boundary_real(:,2),'r')
   
    real_area = polyarea(boundary_real(:,1),boundary_real(:,2));
    disp(['Area of object ' num2str(idx_object) ' = ' num2str(real_area) ' m^2'])


    stats(idx_object).('boundary_pix') = boundary;
    stats(idx_object).('boundary_real') = boundary_real;
    stats(idx_object).('area_real') = real_area;

end 

xlabel('$x \: \rm (m)$')
ylabel('$y \: \rm (m)$')
set(findall(gcf,'-property','FontSize'),'FontSize',15)
set_Papermode(gcf)

figname = [fig_folder prefix_current_fig '_Georeferenced_floes_boundaries_sensitivity' num2str(sensitivity) ...
    '_imopen' num2str(radius_open) '_conn' num2str(conn) '_minarea' num2str(minArea)];
saveas(gcf,[figname '.fig'])
exportgraphics(gcf,[figname '.png'],'Resolution',fig_resolution)
exportgraphics(gcf,[figname '.pdf'],'Resolution',fig_resolution)

%% Organize structure and Save data 
results = struct();
results.stats = stats;
results.BW_param = struct('sensitivity',sensitivity,'radius_open',radius_open,'minArea',minArea,'conn',conn);
results.units = struct('boundary_pix','pix','boundary_real','m','area_real','m');
results.figfilename = img_name;

file2load = [fig_folder prefix_current_fig '_Data_sensitivity' num2str(sensitivity) ...
    '_imopen' num2str(radius_open) '_conn' num2str(conn) '_minarea' num2str(minArea) '.mat'];

disp('Saving Data..')
save(file2load,'results','-v7.3')
disp('Data saved !')






































%% Extract color of ice 

fig = figure; 
imshow(current_img)

% coordinates of buoys
[xi,yi] = getpts(fig);

%%
% select coordinates of a single buoy
x0 = round(xi(1));
y0 = round(yi(1));
w = 150;
crop_img = current_img(y0-w : y0 +w,x0 - w : x0 + w ,:);

% select pixels of buoy
clear xb yb
fig_buoy = figure; 
imshow(crop_img);
[xb,yb] = getpts(fig_buoy);

% P = impixel(im2double(crop_img));

% Average of the colours
idx_x = round(xb);
idx_y = round(yb);

c = zeros(3,length(idx_x));
for i0 = 1:length(idx_x)
    c(:,i0) = squeeze(crop_img(idx_y(i0),idx_x(i0),:)); % colours RGB
end 

mean_RGB = mean(c,2);
disp(mean_RGB)

%% Create an array of distances for each pixel 

% selected values on each channel 
R0 = 138/255; % value on canal R
G0 = 145/255; % value on canal G 
B0 = 165/255; % value on canal B 

img_double = im2double(current_img); % convert image matrix to double
I = sqrt((img_double(:,:,1) - R0).^2 + (img_double(:,:,2) - G0).^2 + (img_double(:,:,3) - B0).^2);
I = ones(size(I)) - I;

figure(1), 
imagesc(I)
colormap("gray")
colorbar()
title('Distance to reference color')

figure(2)
imshow(img_gray)
title('Gray scaled image')

sensitivity = 0.58;
T = adaptthresh(I,sensitivity,'ForegroundPolarity','bright','NeighborhoodSize',2*floor(size(img_gray)/16) + 1);
bw = imbinarize(I,T);
figure(3)
imshow(bw)
title("Binarized image")




%% Contrast enhancement

img_lab = rgb2lab(current_img);
max_luminosity = 100;
L = img_lab(:,:,1)/max_luminosity;


shadow_imadjust = img_lab;
shadow_imadjust(:,:,1) = imadjust(L)*max_luminosity;
shadow_imadjust = lab2rgb(shadow_imadjust);

shadow_histeq = img_lab;
shadow_histeq(:,:,1) = histeq(L)*max_luminosity;
shadow_histeq = lab2rgb(shadow_histeq);

shadow_adapthisteq = img_lab;
shadow_adapthisteq(:,:,1) = adapthisteq(L)*max_luminosity;
shadow_adapthisteq = lab2rgb(shadow_adapthisteq);

figure(1),
imshow(shadow_imadjust)
title('Imadjust')

figure(2),
imshow(shadow_histeq)
title('Histeq')

figure(3),
imshow(shadow_adapthisteq)
title('Adapthisteq')

%%
img_imadjust = imadjust(img_gray);
img_histeq = histeq(img_gray);
img_adapthisteq = adapthisteq(img_gray);

figure(1),
imshow(img_imadjust)
title('Imadjust')

figure(2),
imshow(img_histeq)
title('Histeq')

figure(3),
imshow(img_adapthisteq)
title('Adapthisteq')

%%

% % check Centroid 
% figure, 
% imshow(img_gray)
% hold on 
% plot(new_object.Centroid(1),new_object.Centroid(2),'ro')
% 
% % check BoundingBox
% bbox = [new_object.BoundingBox(2) new_object.BoundingBox(2) + new_object.BoundingBox(4) ; ...
%     new_object.BoundingBox(1) new_object.BoundingBox(1) + new_object.BoundingBox(3)];
% bbox = round(bbox);
% 
% zoom = img_gray(bbox(1,1) : bbox(1,2), bbox(2,1) : bbox(2,2));
% figure
% imshow(zoom)
% 
% % check PixelIdxList
% mask = zeros(sz);
% mask(new_object.PixelIdxList) = 1;
% Ioverlay = labeloverlay(img_gray,mask);
% figure, 
% imshow(Ioverlay)
% 
% % check boundary 
% figure, 
% imshow(img_gray)
% hold on 
% plot(new_object.boundary_pix(:,1),new_object.boundary_pix(:,2),'r.')