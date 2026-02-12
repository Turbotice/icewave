function macro= subimg2main_obj_propreties(micro,topleftxy,sz)
% Convert propreties of a detected object in a sub-image to propreties in
% the main image
% Inupts : - micro_stats : structure, object propreties in the sub-image 
%          - top_leftxy : vector, xy coordinates of the sub-image top left corner in the main image  
%          - sz : size of the main image [ny nx]
% Outputs : - macro_stats : structure containing all propreties of the
% studied object in the main image 

    macro = {};
    
    macro.Centroid = micro.Centroid + topleftxy;
    macro.BoundingBox = zeros([1 4]);
    macro.BoundingBox(1:2) = micro.BoundingBox(1:2) + topleftxy;
    macro.BoundingBox(3:4) = micro.BoundingBox(3:4);
    macro.PixelList = micro.PixelList;
    macro.PixelList(:,1) = round(macro.PixelList(:,1) + topleftxy(1));
    macro.PixelList(:,2) = round(macro.PixelList(:,2) + topleftxy(2));
    macro.PixelIdxList = sub2ind(sz,macro.PixelList(:,2),macro.PixelList(:,1));
    macro.boundary_pix = round(micro.boundary_pix + topleftxy);


end 