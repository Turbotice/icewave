function [results] = ice_floes_detection(img,s,BW_param,fig_folder,prefix_fig,fig_resolution)
% Detect ice floes on a given image and get ice floes area in real space (georeferenced) 
% Inputs :
%   - img : image studied
%   - s : structure containing flight parameters associated to the studied image : 
%       + s.h : drone altitude
%       + s.alpha : angle of the camera to the horizontal 
%       + s.focale : focal length associated to a given format
%       + s.Lx : image length along horizontal (size(img,2))
%       + s.Ly : image length along vertical (size(img,1))
%   - BW_param : structure containing parameters used for binarization and
%   floes detection : 
%       + BW_param.sensitivity : sensitivity used by adaptthresh
%       + BW_param.radius_open : radius of structuring element to perform erosion
%       + BW_param.minArea : minimal area (in pixels) of detected objects
%       + BW_param.conn : connection used for bwboundaries
%   - fig_folder : str, folder where figures will be saved
%   - prefix_fig : str, prefix applied to all figures
%   - fig_resolution : double, resolution with which figures are saved
%
% Outputs :
%   - results : a structure with propreties of detected objects, and associated units, parameters
%   of binarization


% extraction of BW_param :
default_sensitivity = BW_param.sensitivity;
radius_open = BW_param.radius_open;
minArea = BW_param.minArea;
conn = BW_param.conn; 

img_gray = rgb2gray(img);

% define pixel coordinates 
Lx = size(img,2);
Ly = size(img,1);

X_pix=repmat((1:Lx),Ly,1);
Y_pix=repmat((1:Ly)',1,Lx);
% define camera pixel center 
x0 = (Lx + 1)/2;
y0 = (Ly + 1)/2;

% Compute coordinates in real space from pixel coordinates 
[Xreal,Yreal] = projection_real_space(X_pix,Y_pix,x0,y0,s.h,s.alpha,s.focale);

% Apply adaptative threshold to grayscale image
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

    figure,
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

% Open the image
struct_open = strel('disk',radius_open);
BW_cleaned = imopen(BW,struct_open);

% Filter out small objects based on area

BW_filtered = bwareaopen(BW_cleaned, minArea);
figure,
% imshowpair(img_gray,BW_filtered,'montage')
imshow(BW_filtered)
set_Papermode(gcf)
figname = [fig_folder prefix_fig '_BW_filtered_sensitivity' num2str(sensitivity) ... 
    '_imopen' num2str(radius_open) '_conn' num2str(conn) '_minarea' num2str(minArea)];
saveas(gcf,[figname '.fig'])
exportgraphics(gcf,[figname '.png'],'Resolution',fig_resolution)
exportgraphics(gcf,[figname '.pdf'],'Resolution',fig_resolution)

% Detect ROI and extract there position

[B,L,n,~] = bwboundaries(BW_filtered,conn,'noholes','CoordinateOrder',"xy");
disp(['Number of detected objects :' num2str(n)])

stats = regionprops(L,'Centroid','PixelList');

figure,
imshow(img_gray)
hold on 
for k = 1 : n
    plot(stats(k).Centroid(1),stats(k).Centroid(2),'r.')
end 
set_Papermode(gcf)

figname = [fig_folder prefix_fig '_Detected_floes_centroid_sensitivity' num2str(sensitivity) ... 
    '_imopen' num2str(radius_open) '_conn' num2str(conn) '_minarea' num2str(minArea)];
saveas(gcf,[figname '.fig'])
exportgraphics(gcf,[figname '.png'],'Resolution',fig_resolution)
exportgraphics(gcf,[figname '.pdf'],'Resolution',fig_resolution)

% Georeference picture and pixels of a given object + object area computation 
figure,
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

figname = [fig_folder prefix_fig '_Georeferenced_floes_boundaries_sensitivity' num2str(sensitivity) ...
    '_imopen' num2str(radius_open) '_conn' num2str(conn) '_minarea' num2str(minArea)];
saveas(gcf,[figname '.fig'])
exportgraphics(gcf,[figname '.png'],'Resolution',fig_resolution)
exportgraphics(gcf,[figname '.pdf'],'Resolution',fig_resolution)

% Organize structure and Save data 
results = struct();
results.stats = stats;
results.BW_param = struct('sensitivity',sensitivity,'radius_open',radius_open,'minArea',minArea,'conn',conn);
results.units = struct('boundary_pix','pix','boundary_real','m','area_real','m');


end