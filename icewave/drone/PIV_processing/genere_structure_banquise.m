function m = genere_structure_banquise(u,v,u_filt,v_filt,a,filename)

%%%% Change these parameters according to experiment
% fx = 0.1018; %0.0504; % in mm/pixel
% ft = 1/5000; % time scale in s/image
% a : nombre de case à enlever sur les bords, 1
% p & s : parametre de traitements piv
% w : taille de la fenêtre en pixel
%%%%

%w=32; %default value, careful if you change box size in piv processing !

%cel = split(filename,'offset_m');
%cel = split(cel{2},'mV');
%offset = cel{1};
%offset = str2num(offset(1:end));
%m.offset = offset;
%dr(1).z = num2str(str2double(filename(18:20))/50);%relative position in z
%dr(1).h = str2double(filename(11:12)); %liquid height
% 
% disp('Loading data ...')
% load(filename)
% disp('Data loaded')

%%
 
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

%[dimx, dimy]= size(squeeze(u_original{1}));
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
   %if mod(i-2,1000)==0
   %     disp(i-2)
   %end
    u_original{i} = squeeze(u_original{i});
    v_original{i} = squeeze(v_original{i});
    
    m.Vx(:,:,i) = squeeze(u_original{i})';
    m.Vy(:,:,i) = squeeze(v_original{i})';
end

%m.filename = filename;
m.x = (1:1:dimx);
m.y = (1:1:dimy);

dimt = size(m.Vx,3);
m.t = (1:1:dimt);

m.units.Vx = 'pix/frame';
m.units.Vy = 'pix/frame';
m.units.Vz = '';
m.units.x = 'box_idx';
m.units.y = 'box_idx';
m.units.t = 'frame_idx';
     
m.name = filename;%µµstrcat(num2str(dr(1).h), 'mm, z ',num2str(dr(1).z));
m.ysign = +1;
    
%% Remove 2a data points from both x and y 2D plane
[nx,ny,n] = size(m.Vx);

m.Vx = m.Vx(1+a:end-a,1+a:end-a,:);
m.Vy = m.Vy(1+a:end-a,1+a:end-a,:);

m.x = m.x(1:end-2*a);
m.y = m.y(1:end-2*a);

end
