function [Mtab] = sort_points(Mp,nx,ny)
% Sort points of each quadrilateral polygon according to generatecheckerboardpoints ordering
% Inputs : - Mp: structure containing multiple substructures, one for each
% detected polygon. Substructure Mp.m contains field X,Y and corners
%          - nx,ny : ints, number of detected points along each dimension
%          of the polygon
% Output : - Mtab:structure containing all detected points, correctly
% sorted, for each detected polygon / subgrid

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
end