function M = calib_drone_getsubgrid(im,Xout,Yout,name)

figure(1)
hold off
imshow(im)
hold on
plot(Xout,Yout,'ro','markersize',6,'MarkerFaceColor','r')

b=1000;
s=[];
c=0;
n = 9;

M={};

nx = 3840;
ny = 2160;

while isempty(s)
    c=c+1;
    [xc,yc] = getpts(gcf);
    axis([xc-b xc+b,yc-b yc+b])

    plot(xc,yc,'g+','markersize',8,'MarkerFaceColor','g')

    [xc,yc] = getpts(gcf);

    plot(xc,yc,'go','markersize',8,'MarkerFaceColor','g')

    %select axis of the subgrid
    [xh,yh] = getpts(gcf);
    plot(xh,yh,'go','markersize',8,'MarkerFaceColor','g')

    [xv,yv] = getpts(gcf);
    plot(xv,yv,'go','markersize',8,'MarkerFaceColor','g')

    p = [xc,yc];
    u = [xh-xc,yh-yc];
    v = [xv-xc,yv-yc];
    
    eps = 0.1;
    corners{1} = p+(u+v)*(1+1.5*eps);
    corners{2} = p+(u-v)*(1+eps);
    corners{3} = p-(u-v)*(1+1.5*eps);
    corners{4} = p-(u+v)*(1+eps);

    [Xin,Yin] = calib_drone_select_subgrid(Xout,Yout,corners);
    
    plot(Xin,Yin,'bo','markersize',8,'MarkerFaceColor','b')

    M(c).m.Xsub = Xin;
    M(c).m.Ysub = Yin;
    
    save(['transient_' name],'M')
    axis([1,nx,1,ny])
    getframe();

    s=input('Select another subgrid ?')
end
