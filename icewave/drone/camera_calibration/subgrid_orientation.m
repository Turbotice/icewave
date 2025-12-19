%% This scrip enables to correctly oriente subgrids in order to use them for
% camera calibration 

%% Load data 

path2data = ['Z:/skuchly/Codes/' 'matrix_points_bernache.mat'];
load(path2data)

%% Show a given subgrid
nx = 9;
ny = 7;

for i = 1:length(Mp)
    disp(i)
    X = Mp(i).m.X;
    Y = Mp(i).m.Y;
    
    points = [X(:),Y(:)];
    K = convhull(points(:,1),points(:,2)); % convex hull
    % Hull points
    Hullpts = points(K(1:end-1,:),:);

% reduce number of hull points using interior angle of each hull point. 

    corners = Hullpts;
    
    while size(corners,1) > 4
        N = size(corners,1);
        angles = zeros(N,1);
    
        for k = 1:N
            p_prev = corners(mod(k-2,N)+1,:);
            p_curr = corners(k,:);
            p_next = corners(mod(k,N)+1,:);
    
            v1 = p_prev - p_curr;
            v2 = p_next - p_curr;
    
            % Interior angle
            angles(k) = acos(dot(v1,v2)/(norm(v1)*norm(v2)));
        end
    
        % Remove point with largest angle (least corner-like)
        [~, idx] = max(angles);
        corners(idx,:) = [];
    end


    % Plot subgrid and detected points 
    
    % figure(1)
    % plot(points(:,1),points(:,2),'+')
    % hold on
    % plot(Hullpts(:,1),Hullpts(:,2),'o')
    % hold on 
    % plot(corners(:,1),corners(:,2),'ks')
    % text(corners(:,1)+2,corners(:,2)+2,{'1','2','3','4'},'Color','k')
    
    
    % Define ordering of corners
    
    % corners_order = {1,2,3,4;'top_left','top_right','bottom-right','bottom-left'};
    
    % compute centroid coordinates
    C = mean(corners,1);
    
    % compute angle to centroid
    angles = atan2(corners(:,2)- C(2),corners(:,1) - C(1));
    
    % sort counterclockwise
    [~,idx] = sort(angles);
    corners_ccw = corners(idx,:);
    
    % enforce starting point as top left
    % compute angle to centroid
    angles_ccw = atan2(corners_ccw(:,2)- C(2),corners_ccw(:,1) - C(1));
    % [~,topleft] = min(corners_ccw(:,1) + corners_ccw(:,2));
    [~,topleft] = max(angles_ccw);
    ordered = circshift(corners_ccw,1-topleft);
    % invert corner 2 and 4
    ordered([2,4],:) = ordered([4,2],:);

    % Find points between two corners 
    T = 10;
    % line between corners 1 and 2
    [x_line,y_line] = find_line(X,Y,ordered(1,1),ordered(1,2),...
        ordered(2,1),ordered(2,2),T);
    
    nb_points = length(x_line);
    if nb_points == ny
        ordered([2,4],:) = ordered([4,2],:);
    end

    % Visualisation 
    
    fig = figure(2);
    plot(points(:,1),points(:,2),'+')
    hold on 
    plot(ordered(:,1),ordered(:,2),'ro')
    text(ordered(:,1)+2,ordered(:,2)+2,{'1','2','3','4'},'Color','k')
    hold off
    % plot(x_line,y_line,'diamond','Color','green')
    title(i)

    pause(0.1)
    % plot(corners_ccw(:,1),corners_ccw(:,2),'ro')
    % text(corners_ccw(:,1)+2,corners_ccw(:,2)+2,{'1','2','3','4'},'Color','k')

    % Store corners 
    Mp(i).m.corners = ordered;

end 

%% Sort all points of a subgrid according to generatecheckerboardpoints ordering

Mtab = {};
c = 1;
T = 10;
for i = 1:length(Mp)
    disp(i)
    corners = Mp(i).m.corners;
    X = Mp(i).m.X;
    Y = Mp(i).m.Y;
    points = [X(:),Y(:)];
    ordered_pts = points;
    
    % find line between corner 1 and 2
    [xtop,ytop] = find_line(X,Y,corners(1,1),corners(1,2),...
        corners(2,1),corners(2,2),T);
    
    % find line between corner 4 and 3
    [xbottom,ybottom] = find_line(X,Y,corners(4,1),corners(4,2),...
        corners(3,1),corners(3,2),T);

    if and(length(xtop) == nx,length(xbottom) == nx)
        
        % sort each column of the grid 
        for k = 1:nx
            c_top = [xtop(k),ytop(k)];
            c_bott = [xbottom(k),ybottom(k)];
            [xline,yline] = find_line(X,Y,xtop(k),ytop(k),...
                xbottom(k),ybottom(k),T);
            
            ordered_pts(1 + ny*(k-1):ny*k,1) = xline;
            ordered_pts(1 + ny*(k-1):ny*k,2) = yline;
        
            figure(4),
            plot(points(:,1),points(:,2),'+')
            hold on 
            plot(xtop,ytop,'ro')
            hold on 
            plot(xbottom,ybottom,'dg')
            hold on 
            plot(xline,yline,'ks')
            title(i)
            hold off
        
            % pause(0.1)
        end
    
        figure(5),
        plot(ordered_pts(:,1),ordered_pts(:,2),'+') 
        for k = 1 : size(ordered_pts,1)
            text(ordered_pts(k,1)+2,ordered_pts(k,2)+2,num2str(k))
        end
    
        Mtab(c).points = ordered_pts;
        c = c + 1;
    end
end 

%% Save Mtab

path2save = 'Y:/Banquise/Sebastien/Calibration_intrinseque/Calib_PMMH/';
filename = [path2save 'matrix_image_points_bernache.mat'];

save(filename,'Mtab','-v7.3')

%% Estimate camera calibration points 

% collect imagepoints
for i = 1:length(Mtab)
    imagePoints(:,:,i) = Mtab(i).points;
end

% generate checkerboard pattern 
boardSize = [8,10];
squareSizeInMM = 3955 / 40;
worldPoints = generateCheckerboardPoints(boardSize,squareSizeInMM);
imageSize = [3840,2160];

% calibration
stereoParams = estimateCameraParameters(imagePoints,worldPoints,ImageSize=imageSize);



%%

i = 31;
corners = Mp(i).m.corners;
X = Mp(i).m.X;
Y = Mp(i).m.Y;
points = [X(:),Y(:)];
ordered_pts = points;
T = 10;

% find line between corner 1 and 2
[xtop,ytop] = find_line(X,Y,corners(1,1),corners(1,2),...
    corners(2,1),corners(2,2),T);

% find line between corner 4 and 3
[xbottom,ybottom] = find_line(X,Y,corners(4,1),corners(4,2),...
    corners(3,1),corners(3,2),T);


figure(4),
plot(points(:,1),points(:,2),'+')
hold on 
plot(corners(:,1),corners(:,2),'ks')
hold on 
plot(xtop,ytop,'ro')
hold on 
plot(xbottom,ybottom,'dg')









%% Possible solution system coordinate sorting

% define vectors
u = corners(2,:) - corners(1,:);   % horizontal direction
v = corners(4,:) - corners(1,:);   % vertical direction

u = u / norm(u);
v = v / norm(v);

% project all points
P0 = points - corners(1,:);   % shift origin to top-left corner

s = P0 * u';   % horizontal coordinate
t = P0 * v';   % vertical coordinate

% normalize coordinate to grid indices
s = (s - min(s)) / (max(s) - min(s));
t = (t - min(t)) / (max(t) - min(t));

col = round(s * (nx-1));
row = round(t * (ny-1));