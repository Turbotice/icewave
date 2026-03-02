
function [Xin,Yin] = calib_drone_select_subgrid(X,Y,corners)


for j=1:length(corners)
    plot(corners{j}(1),corners{j}(2),'ks','markersize',10,'MarkerFaceColor','k')
end
getframe();

% click on image to select points 

% select only points in the rectangl defined by tape
xv = [];
yv = [];
for j = 1 : length(corners)
    xv = cat(2,xv,corners{j}(1));
    yv = cat(2,yv,corners{j}(2));
end 
xv([3 4]) = xv([4 3]);
yv([3 4]) = yv([4 3]);
in = inpolygon(X,Y,xv,yv);

%[vx,vy] = voronoi(X(in),Y(in));
%voronoi(X(in),Y(in))

Xin = X(in);
Yin = Y(in);