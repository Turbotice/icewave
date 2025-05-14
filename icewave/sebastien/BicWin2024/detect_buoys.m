function new_table = detect_buoys(img,RGB0,threshold,SE,min_distance)
    % Detect object on an image based on their colour
    % Inputs : - img, image to work on 
    %          - RGB0, tuple of reference colour. Contains RGB values of
    %          objects we want to detect
    %          - threshold, value (between 0 and 1) of threshold above which a pixel is
    %          detected
    %          - SE, matlab structuring element to perform erosion
    %          - min_distance, integer, minimal distance (in pixel) between pixels that
    %          need to be merged
    % Outputs : - new_table : table containing all properties of detected
    % objects 

    % Build a distance from colours
    img_double = im2double(img);
    % define center for color coordinates system
    R0 = RGB0(1);
    G0 = RGB0(2);
    B0 = RGB0(3);
    % defin color intensity in color coordinates system
    I = sqrt((img_double(:,:,1) - R0).^2 + (img_double(:,:,2) - G0).^2 + (img_double(:,:,3) - B0).^2);
    % Keep pixels corresponding to the indicated color
    Inew = ones(size(I)) - I;

    Inew(Inew >= threshold) = 1;
    Inew(Inew < threshold) = 0;

    % Detect ROI and extract there position 

    Inew = imerode(Inew,SE); % erosion 
    CC = bwconncomp(Inew);

    stats = regionprops("table",CC);
    % Plot detected objects
    stats = sortrows(stats,'Centroid','descend'); % sort rows according to x position
  
    % Compute distance matrix between objects 
    D = squareform(pdist(stats.Centroid)); 
    
    % Merge particles that are too close to each other 
    mask_dist = (D < min_distance) & (D > 0);
    
    if ismember(1,mask_dist)
        [r,c] = find(mask_dist); % find indices of non-zero elements
        % select only superior part of the matrix 
        list_i = r(r>c);
        list_j = c(r>c);

        particle_i = list_i(~ ismember(list_i,list_j));
        particle_j = list_j(~ ismember(list_i,list_j));

        unique_part = unique(particle_i);

        new_part_centroid = zeros(length(unique_part),2);
        new_part_area = zeros(length(unique_part),1);
        for i0 = 1:length(unique_part)
            part = unique_part(i0);
            associated_part = particle_j(particle_i == part);
            list_particle = cat(1,part,associated_part);

            % computes mass center coordinates 
            new_part_centroid(i0,:) = sum(stats.Area(list_particle).*stats.Centroid(list_particle,:))...
                /sum(stats.Area(list_particle));
            new_part_area(i0) = sum(stats.Area(list_particle));
        end 
        % Get particles that do not need to be merged

        idx_list = (1:size(D,1));
        non_merged_mask = ~ (ismember(idx_list,particle_i) + ismember(idx_list,particle_j));
        non_merged = idx_list(non_merged_mask);

        non_merged_centroid = stats.Centroid(non_merged,:);
        non_merged_area = stats.Area(non_merged);
        new_centroid = cat(1,non_merged_centroid,new_part_centroid);
        new_area = cat(1,non_merged_area,new_part_area);
        
        Centroid = new_centroid;
        Area = new_area;
        new_table = table(Centroid,Area);
    else 
        
        new_table = stats;
    end 
    
    new_table = sortrows(new_table,'Centroid','descend');
    
end 