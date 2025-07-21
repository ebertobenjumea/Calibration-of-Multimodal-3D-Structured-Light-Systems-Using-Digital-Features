close all
clc
clear all
direc='D:\R3D\Adquisiciones\Experimento 54\Objects\SL\';
figure;
error_rms=[];
for i=8:28
    if i>9
        aux=num2str(i);
    else
        aux=['0' num2str(i)];
    end
    
% Cargamos digital features image
    pc_dig = pcread([direc 'Digital Features 1 UTB\new_obj_vis_ir_' aux '.ply']);

%     subplot(121)
%     pcshow(pc_dig,'BackgroundColor','w');
%     xlabel('x (mm)')
%     ylabel('y (mm)')
%     zlabel('z (mm)')
%     c = colorbar; 
%     c.Label.String = '° C';
%     set(gca, 'Zdir', 'reverse')
%     set(gca, 'Xdir', 'reverse')
%     view(65,42)
%     grid off
%     %zlim([350 430])
%     title('Digital features I')
    
    % Cargamos traditional image
    pc_trad = pcread([direc 'Two Stages 2\new_obj_vis_ir_' aux '.ply']);
%     subplot(122)
%     pcshow(pc_trad,'BackgroundColor','w');
%     xlabel('x (mm)')
%     ylabel('y (mm)')
%     zlabel('z (mm)')
%     c = colorbar; 
%     c.Label.String = '° C';
%     set(gca, 'Zdir', 'reverse')
%     set(gca, 'Xdir', 'reverse')
%     view(65,42)
%     
%     grid off
%     %zlim([350 430])
%     title('Traditional')
%     sgtitle(['Pose ' aux])

    % Creamos matrices de información
    % Digital features I
    x_dig=pc_dig.Location(:,1);
    y_dig=pc_dig.Location(:,2);
    z_dig=pc_dig.Location(:,3);
    t_dig=pc_dig.Intensity;
    %[X_dig,Y_dig]=meshgrid([min(x_dig):0.05:max(x_dig)], [min(y_dig):0.05:max(y_dig)]);
    
    [X_dig,Y_dig]=meshgrid(linspace(min(x_dig),max(x_dig),1280), linspace(min(y_dig),max(y_dig),1024));
    
    
    
    z_est_dig=griddata(x_dig,y_dig,z_dig,X_dig,Y_dig,'cubic');
    t_est_dig=griddata(x_dig,y_dig,t_dig,X_dig,Y_dig,'cubic');
    z_est_dig(z_est_dig>max(z_dig))=NaN;
    t_est_dig(t_est_dig>max(t_dig))=NaN;
    z_est_dig(z_est_dig<min(z_dig))=NaN;
    t_est_dig(t_est_dig<min(t_dig))=NaN;
    
    % Visualización de imagenes 3D y térmicas
    
    
    % Visualizacion 3D
    % figure;
    subplot(121)
    s = surf(-X_dig,Y_dig,-z_est_dig,'FaceColor', 'interp',...
                'EdgeColor', 'none',...
                'FaceLighting', 'phong');
    
    
    set(gca, 'DataAspectRatio', [1, 1, 1])
    axis equal;
    view(65,42)
    camlight right
    axis on
    grid on
    xlabel('x (mm)')
    ylabel('y (mm)')
    zlabel('z (mm)')
    
    
    
    
    disp('Cantidad de datos NaN:')
    sum(sum(isnan(z_est_dig))) % Ojo con esto
    sum(sum(isnan(t_est_dig))) % Ojo con esto
    
    
    
    % Traditional
    x_trad=pc_trad.Location(:,1);
    y_trad=pc_trad.Location(:,2);
    z_trad=pc_trad.Location(:,3);
    t_trad=pc_trad.Intensity;
    %[X_trad,Y_trad]=meshgrid([min(x_trad):0.05:max(x_trad)], [min(y_trad):0.05:max(y_trad)]);
    [X_trad,Y_trad]=meshgrid(linspace(min(x_trad),max(x_trad),1280), linspace(min(y_trad),max(y_trad),1024));
    
    
    
    z_est_trad=griddata(x_trad,y_trad,z_trad,X_trad,Y_trad,'cubic');
    t_est_trad=griddata(x_trad,y_trad,t_trad,X_trad,Y_trad,'cubic');
    z_est_trad(z_est_trad>max(z_trad))=NaN;
    t_est_trad(t_est_trad>max(t_trad))=NaN;
    z_est_trad(z_est_trad<min(z_trad))=NaN;
    t_est_trad(t_est_trad<min(t_trad))=NaN;
    
    
    
    % Visualizacion 3D
    subplot(122);
    
    
    s = surf(-X_trad,Y_trad,-z_est_trad,'FaceColor', 'interp',...
                 'EdgeColor', 'none',...
                 'FaceLighting', 'phong');
    set(gca, 'DataAspectRatio', [1, 1, 1])
    axis equal;
    view(65,42)
    camlight right
    axis on
    grid on
    xlabel('x (mm)')
    ylabel('y (mm)')
    zlabel('z (mm)')

    
    disp('Cantidad de datos NaN:')
    sgtitle(['Pose ' aux])
%    Error
    %tamano=300:700;
    %X_dig=X_dig(tamano,tamano);
    %Y_dig=Y_dig(tamano,tamano);
    %z_est_dig=z_est_dig(tamano,tamano);
    x=X_dig(:);
    y=Y_dig(:);
    z=z_est_dig(:);
    map_z=~isnan(z);
    x_n=x(map_z);
    y_n=y(map_z);
    z_n=z(map_z);
%     sum(sum(isnan(x_n)))
%     sum(sum(isnan(y_n)))
%     sum(sum(isnan(z_n)))

    DM = [x_n, y_n, ones(size(z_n))];                             % Design Matrix
                                            % Estimate Parameters
    %B=pinv(DM,z_n);
    B=lsqminnorm(DM,z_n,'warn');
    z_hat = B(1)*x_n + B(2)*y_n + B(3)*ones(size(x_n));
    error=z_hat-z_n;
    z_est=double(isnan(z));
    z_est(z_est==1)=NaN;
    %z_est=z;
    error_2d=z;
    cnt_data=0;
    m_mask=zeros(size(z));
    for j=1:length(z_est)
        if ~isnan(z_est(j))
            cnt_data=cnt_data+1;
            z_est(j)=z_hat(cnt_data);
            error_2d(j)=error(cnt_data);
            m_mask(j)=1;
        end
    end
    %map_error=reshape(error_2d,[length(tamano),length(tamano)]);
    %z_est=reshape(z_est,[length(tamano),length(tamano)]);
        map_error=reshape(error_2d,[1024,1280]);
    z_est=reshape(z_est,[1024,1280]);

    subplot(232)
    s = surf(X_dig,Y_dig,z_est_dig,'FaceColor', 'interp',...
                    'EdgeColor', 'none',...
                    'FaceLighting', 'phong');

    set(gca, 'DataAspectRatio', [1, 1, 1])
    axis equal;
    view(30,30);
    camlight right
    axis on
    grid on
    xlabel('x')
    ylabel('y')
    zlabel('z')
    title('3D ideal plane')

    subplot(235)
    s = surf(map_error,'FaceColor', 'interp',...
                    'EdgeColor', 'none',...
                    'FaceLighting', 'phong');

    subplot(236)
    imagesc(map_error)
    colorbar
    title('Z error')

    subplot(234)
    %m_mask=reshape(m_mask,[length(tamano),length(tamano)]);
        m_mask=reshape(m_mask,[1024,1280]);
    imshow(m_mask)
    histogram(error)
    title('Hist')
    %imhist(map_error)
    sgtitle(['Pose ' num2str(i-1)])
    error_rms(i)=sqrt((sum(error.^2))/length(error));
    disp(['RMSE: ' num2str(error_rms(i))])

    pause(1)

end

figure;
plot(error_rms([8:20 22:end]))
xlabel('Pose')
ylabel('RMSE')
title('RMSE')

mean_error=mean(error_rms([8:20 22:end]));
desv=std(error_rms([8:20 22:end]));
disp(['RMSE= ' num2str(mean_error) '±' num2str(desv)])

figure; 
    s = surf(map_error,'FaceColor', 'interp',...
                    'EdgeColor', 'none',...
                    'FaceLighting', 'phong');

%%

close all
clc
figure();

%n_poses=50;
error_rms=[];
for i=inicio:n_poses
    disp(['Pose ' num2str(i-1)])
    if i<11
        aux=['0' num2str(i-1)];
    else
        aux=num2str(i-1);
    end
    
    load([direc camera1 'r3d_smooth_pose_' aux '.mat'])
    width=size(ZcM,2);
    height=size(ZcM,1);
    subplot(231)
    s = surf(XcM,YcM,ZcM,'FaceColor', 'interp',...
                    'EdgeColor', 'none',...
                    'FaceLighting', 'phong');

    set(gca, 'DataAspectRatio', [1, 1, 1])
    axis equal;
    view(30,30);
    camlight right
    axis on
    grid on
    xlabel('x')
    ylabel('y')
    zlabel('z')
    title('3D data')
    
    subplot(233)
    imagesc(ZcM)
    colorbar
    title('Z')
% 
%     sum(sum(isnan(XcM)))
%     sum(sum(isnan(YcM)))
%     sum(sum(isnan(ZcM)))
%     x=XcM(:);
%     y=YcM(:);
    [x,y]=meshgrid(1:width,1:height);
    x=x(:);
    y=y(:);
    z=ZcM(:);
    map_z=~isnan(z);
    x_n=x(map_z);
    y_n=y(map_z);
    z_n=z(map_z);
%     sum(sum(isnan(x_n)))
%     sum(sum(isnan(y_n)))
%     sum(sum(isnan(z_n)))

    DM = [x_n, y_n, ones(size(z_n))];                             % Design Matrix
                                            % Estimate Parameters
    %B=pinv(DM,z_n);
    B=lsqminnorm(DM,z_n,'warn');
    z_hat = B(1)*x_n + B(2)*y_n + B(3)*ones(size(x_n));
    error=z_hat-z_n;
    %map_error=reshape(error,[height,width]);

%      subplot(235)
%      plot(error,'.')
%      title('Error')




%     subplot(223)
%     plot(z_n,'.')
%     subplot(224)
%     plot(z_hat,'.')
    z_est=double(isnan(z));
    z_est(z_est==1)=NaN;
    %z_est=z;
    error_2d=z;
    cnt_data=0;
    m_mask=zeros(size(z));
    for j=1:length(z_est)
        if ~isnan(z_est(j))
            cnt_data=cnt_data+1;
            z_est(j)=z_hat(cnt_data);
            error_2d(j)=error(cnt_data);
            m_mask(j)=1;
        end
    end
    map_error=reshape(error_2d,[height,width]);
    z_est=reshape(z_est,[height,width]);

    subplot(232)
    s = surf(XcM,YcM,z_est,'FaceColor', 'interp',...
                    'EdgeColor', 'none',...
                    'FaceLighting', 'phong');

    set(gca, 'DataAspectRatio', [1, 1, 1])
    axis equal;
    view(30,30);
    camlight right
    axis on
    grid on
    xlabel('x')
    ylabel('y')
    zlabel('z')
    title('3D ideal plane')

    subplot(235)
    s = surf(map_error,'FaceColor', 'interp',...
                    'EdgeColor', 'none',...
                    'FaceLighting', 'phong');

    subplot(236)
    imagesc(map_error)
    colorbar
    title('Z error')

    subplot(234)
    m_mask=reshape(m_mask,[height,width]);
    imshow(m_mask)
    histogram(error)
    title('Hist')
    %imhist(map_error)
    sgtitle(['Pose ' num2str(i-1)])
    error_rms(i)=sqrt((sum(error.^2))/length(error));
    disp(['RMSE: ' num2str(error_rms(i))])
    pause(1);

end
figure;
plot(error_rms)
xlabel('Pose')
ylabel('RMSE')
title('RMSE')

mean_error=mean(error_rms(inicio:n_poses));
desv=std(error_rms(inicio:n_poses));
disp(['RMSE= ' num2str(mean_error) '±' num2str(desv)])

figure; 
    s = surf(map_error,'FaceColor', 'interp',...
                    'EdgeColor', 'none',...
                    'FaceLighting', 'phong');
figure; 
    s = surf(Fx,'FaceColor', 'interp',...
                    'EdgeColor', 'none',...
                    'FaceLighting', 'phong');


    figure;
    plot(Fy(600,400:1000))


