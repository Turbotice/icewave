function [M] = find_corners(Mp,nx,ny)
% Find corners pixel coordinates of a convex quadrilateral polygon
% Inputs: - Mp, structure containing several substructures, one for
% each detected polygon. Mp.m is a
% substructure with fields Mp.m.X and Mp.m.Y which are coordinates of
% detected tiles / squares
% size of Mp.m.X and Mp.m.Y is nx*ny
%         - nx,ny = ints, number of points along each dimensions of the
%         quadrilateral polygon

    for i = 1:length(Mp) % loop over all detected polygons
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
    
    M = Mp;

end