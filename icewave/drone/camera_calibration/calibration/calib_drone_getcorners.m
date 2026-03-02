function corners =calib_drone_getcorners(im)

figure(1)
hold off
imshow(im)
for i=1:4
    [x,y] = getpts(gcf);
    regions{i} = [x,y];
end

b = 100;
for i=1:4
    xc=regions{i}(1);
    yc=regions{i}(2);
    
    axis([xc-b xc+b,yc-b yc+b])

   [x,y] = getpts(gcf);
   corners{i} = [x,y];
end
