function plot_velocity_features(Vx,Vy,fx,i,caxis_amp,fig_folder)

% This function plots average and standard deviation of velocity fields 

% Arguments : 
% - Vx:  velocity field along x obtained from PIV process
% - Vy : velocity field along x obtained from PIV process
% - fx : space acquisition frequency 
% - i : index of the frame shown
% - caxis_amp : amplitude of the colorbar
% - fig_folder : folder where to save plots 

    [nx,ny,nt] = size(Vx);

    % create a meshgrid
    y = (ny:-1:1);
    x = (1:1:nx);
    [X,Y]=meshgrid(x,y);

    % % compute the mean component of the velocity field for each frame
    % Vxmoy = mean(mean(m.Vx,2),1);
    % Vymoy = mean(mean(m.Vy,2),1);
    % % reduce the velocity field from its mean value
    % m.Vx = m.Vx - Vxmoy;
    % m.Vy = m.Vy - Vymoy;

    %% get an image of the velocity fields
    figure;
    
    subplot(2,1,1)
    hold off
    surf(X./fx,Y./fx,Vx(:,:,i)')
    shading interp
    view(2)
    axis([0 size(X,2)/fx 0 size(X,1)/fx]);
    caxis([-caxis_amp caxis_amp])
    cbar =colorbar();
    cbar.Label.String = '$V_x \: \rm (m.s^{-1})$';
    cbar.Label.Interpreter = 'latex';
    cbar.FontSize = 13;
    %title('$V_x$','Interpreter','latex')
    xlabel('$x$ (m)','Interpreter','latex');
    ylabel('$y$ (m)','Interpreter','latex');
    ax = gca;
    ax.FontSize = 13;

    subplot(2,1,2)
    hold off
    surf(X./fx,Y./fx,Vy(:,:,i)')
    axis([0 size(X,2)/fx 0 size(X,1)/fx]);
    shading interp
    view(2)
    caxis([-caxis_amp caxis_amp])
    cbar =colorbar();
    cbar.Label.String = '$V_y \: \rm (m.s^{-1})$';
    cbar.Label.Interpreter = 'latex';
    cbar.FontSize = 13;
    %title('$V_y$','Interpreter','latex')
    xlabel('$x$ (m)','Interpreter','latex');
    ylabel('$y$ (m)','Interpreter','latex');
    ax = gca;
    ax.FontSize = 13;

    % set correctly the image position for a pdf format 
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(gcf,'filename','-dpdf','-r0')

    fig_file = [fig_folder 'Velocity_map_idx_' num2str(i) 'caxis_amp_' num2str(caxis_amp)];
    saveas(gcf,fig_file,'fig');

    %% Moyenne, ecart type

    % computes mean and standard deviation over time 
    Vxmoy = mean(Vx,3,'omitnan');
    Vymoy = mean(Vy,3,'omitnan');
    Vxstd = std(Vx,[],3);
    Vystd = std(Vy,[],3);

    fig_average = figure(2);
    subplot(2,1,1)
    surf(X./fx,Y./fx,Vxmoy')
    %title('$\langle V_x \rangle _t$','Interpreter','latex','FontSize',13)
    xlabel('$x$ (m)','Interpreter','latex');
    ylabel('$y$ (m)','Interpreter','latex');
    shading interp
    view(2)
    axis([0 size(X,2)/fx 0 size(X,1)/fx]);
    % caxis([-3 3])
    cbar = colorbar();
    cbar.Label.String = '$\langle V_x \rangle _t \: \rm (m.s^{-1})$';
    cbar.Label.Interpreter = 'latex';
    cbar.FontSize = 13;
    ax = gca;
    ax.FontSize = 13;

    subplot(2,1,2)
    surf(X./fx,Y./fx,Vymoy')
    %title('$\langle V_y \rangle _t$','Interpreter','latex')
    xlabel('$x$ (m)','Interpreter','latex');
    ylabel('$y$ (m)','Interpreter','latex');
    shading interp
    view(2)
    axis([0 size(X,2)/fx 0 size(X,1)/fx]);
    % caxis([-2 2].*scale_V)
    cbar = colorbar();
    cbar.Label.String = '$\langle V_y \rangle _t \: \rm (m.s^{-1})$';
    cbar.Label.Interpreter = 'latex';
    cbar.FontSize = 13;
    ax = gca;
    ax.FontSize = 13;    

    % set correctly the image position for a pdf format 
    set(fig_average,'Units','Inches');
    pos = get(fig_average,'Position');
    set(fig_average,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(fig_average,'filename','-dpdf','-r0')

    fig_std = figure(3);
    subplot(2,1,1)
    surf(X./fx,Y./fx,Vxstd')
    %title('$\sigma (V_x)$','Interpreter','latex')
    xlabel('$x$ (m)','Interpreter','latex');
    ylabel('$y$ (m)','Interpreter','latex');
    shading interp
    view(2)
    axis([0 size(X,2)/fx 0 size(X,1)/fx]);
    % caxis([-3 3])
    cbar = colorbar();
    cbar.Label.String = '$\sigma (V_x) \: \rm (m.s^{-1})$';
    cbar.Label.Interpreter = 'latex';
    cbar.FontSize = 13;
    ax = gca;
    ax.FontSize = 13;

    subplot(2,1,2)
    surf(X./fx,Y./fx,Vystd')
    %title('$\sigma (V_y)$','Interpreter','latex')
    xlabel('$x$ (m)','Interpreter','latex');
    ylabel('$y$ (m)','Interpreter','latex');
    shading interp
    view(2)
    axis([0 size(X,2)/fx 0 size(X,1)/fx]);
    % caxis([-3 3])
    cbar = colorbar();
    cbar.Label.String = '$\sigma (V_y) \: \rm (m.s^{-1})$';
    cbar.Label.Interpreter = 'latex';
    cbar.FontSize = 13;
    ax = gca;
    ax.FontSize = 13;

    % set correctly the image position for a pdf format 
    set(fig_std,'Units','Inches');
    pos = get(fig_std,'Position');
    set(fig_std,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(fig_std,'filename','-dpdf','-r0')

    fig_file_mean = [fig_folder 'Velocity_average'];
    fig_file_std = [fig_folder 'Velocity_std'];

    saveas(fig_average, fig_file_mean, 'fig');
    saveas(fig_std, fig_file_std, 'fig');

end 
