function [Xlist,Ylist]=calib_drone_correcpoints(m)

b=50;

corners = m.corners;
xv = [];
yv = [];
for j = 1 : length(corners)
    xv = cat(2,xv,corners{j}(1));
    yv = cat(2,yv,corners{j}(2));
end 

xmin = min(xv);
xmax = max(xv);
ymin = min(yv);
ymax = max(yv);

nx = 3840;
ny = 2160;

s=[];
c=0;
while isempty(s)
    c=c+1;
    [xc,yc] = getpts(gcf);
    axis([xc-b xc+b,yc-b yc+b])

    [x,y] = getpts(gcf);

    Xlist(c)=x;
    Ylist(c)=y;

    plot(x,y,'bo','markersize',8,'MarkerFaceColor','b')

    axis([1,nx,1,ny])
    s=input('Look for other points ?')
end

