function m = genere_structure_banquise(u,v,u_filt,v_filt,xpix,ypix,a,filename)

m={};
if exist('u')
    u_original = squeeze(u);
    v_original = squeeze(v);
end
if exist('u_filt')
    if ~isempty(u_filt{1})
        disp('use filtered data')
        u_original = u_filt;
        v_original = v_filt;
    else
        if exist('u_original')
            u_original = squeeze(u_original);
            v_original = squeeze(v_original);
        end
    end
end

c = length(u_original);
for i=1:c
    [dimx, dimy]= size(squeeze(u_original{i}));
    if dimx>0
        n=i;
    end
end

[dimy, dimx]= size(squeeze(u_original{1}));
m.Vx = zeros(dimx,dimy,n); % Remove 2b data points from both x and y
m.Vy = zeros(dimx,dimy,n); % Remove 2b data points from both x and y
m.Vz = ''; % create an empty field 

for i=1:n

    u_original{i} = squeeze(u_original{i});
    v_original{i} = squeeze(v_original{i});
    
    m.Vx(:,:,i) = squeeze(u_original{i})';
    m.Vy(:,:,i) = squeeze(v_original{i})';
end

m.x = (1:1:dimx);
m.y = (1:1:dimy);
[X,Y] = meshgrid(m.x,m.y);
m.X = permute(X,[2,1]);
m.Y = permute(Y,[2,1]);

dimt = size(m.Vx,3);
m.t = (1:1:dimt);

m.units.Vx = 'pix/frame';
m.units.Vy = 'pix/frame';
m.units.Vz = '';
m.units.x = 'box_idx';
m.units.y = 'box_idx';
m.units.X = 'box_idx';
m.units.Y = 'box_idx';
m.units.t = 'frame_idx';
     
m.name = filename;

% create a substructure for boxes PIXEL coordinates 
m.PIXEL = struct('x_pix',xpix,'y_pix',ypix);
    
%% Remove 2a data points from both x and y 2D plane
[nx,ny,n] = size(m.Vx);

m.Vx = m.Vx(1+a:end-a,1+a:end-a,:);
m.Vy = m.Vy(1+a:end-a,1+a:end-a,:);

m.x = m.x(1:end-2*a);
m.y = m.y(1:end-2*a);

m.PIXEL.x_pix = m.PIXEL.x_pix(1 + a : end - a);
m.PIXEL.y_pix = m.PIXEL.y_pix(1 + a : end - a);


end
