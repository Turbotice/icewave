function obj_ind = find_objects(x,y,stats,sz)
    % Find objects within a regionprops structure corresponding to a set of
    % (x,y) coordinates
    % Inputs : - x : vector of x coordinates (horizontal of image)
    %          - y : vector of y coordinates (vertical of image)
    %          - stats : structure of pixels belonging to each objects
    %          (usually extracted from regionprops), it should have a field
    %          : 'PixelIdxList'
    %          - sz : size of the 2D image
    % Output : - obj_ind : indices of objects found

    head = 'PixelIdxList';
    if ~isfield(stats,head)
        error('Structure of objects has no field called PixelIdxList')
    end 
    
    
    x = floor(x);
    y = floor(y);
    linear_idx = sub2ind(sz,y,x);
    obj_ind = zeros(size(x));
    
    for i = 1 : length(x)
        for j = 1:length(stats)
            if ismember(linear_idx(i),stats(j).PixelIdxList)
                obj_ind(i) = j;
            end 
        end
    end

end 