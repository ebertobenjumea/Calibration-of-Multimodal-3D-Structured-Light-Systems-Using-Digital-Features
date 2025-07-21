close all
clc
clear all
direc='D:\R3D\Adquisiciones\Experimento 54\Objects\SL\';
figure;
error_rms=[];
%for i=11:11
    i=11
    if i>9
        aux=num2str(i);
    else
        aux=['0' num2str(i)];
    end
    
% Cargamos digital features image
    pc_dig = pcread([direc 'Digital Features 1 UTB\new_obj_vis_ir_' aux '.ply']);
    figure; 
    vis_texture=imread([direc 'pose_' aux '\im_51.png']);
    imshow(vis_texture)
    level = graythresh(vis_texture);        %Umbralizamos
    BW=imbinarize(vis_texture,level);
%     figure;
%     imshow(BW)
    mirror_corners = find_corners_mirror(BW,1280);
    mirror_corners(1,1)=mirror_corners(1,1)+10;
    mirror_corners(1,4)=mirror_corners(1,4)+10;
    mirror_corners(1,2)=mirror_corners(1,2)-10;
    mirror_corners(1,3)=mirror_corners(1,3)-10;
    mirror_corners(2,1)=mirror_corners(2,1)+10;
    mirror_corners(2,4)=mirror_corners(2,4)-10;
    mirror_corners(2,2)=mirror_corners(2,2)+10;
    mirror_corners(2,3)=mirror_corners(2,3)-10;
    hold on
    %plot(mirror_corners(1,:),mirror_corners(2,:),'+r')
    line([mirror_corners(1,:) mirror_corners(1,1)],[mirror_corners(2,:) mirror_corners(2,1)],'Color','blue','LineStyle','--','Linewidth',1.5)
    %BW = zeros([1024,1280]);
    mirror_corners=double(mirror_corners);
    
    % Usar la función poly2mask para crear la máscara a partir de las coordenadas
    BW = poly2mask(mirror_corners(1,:)', mirror_corners(2,:)', 1024, 1280);
    figure;
    imshow(BW)
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
    % Creamos matrices de información
    % Digital features I
    x_dig=pc_dig.Location(:,1);
    y_dig=pc_dig.Location(:,2);
    z_dig=pc_dig.Location(:,3);
    t_dig=pc_dig.Intensity;
    mask_dig=load([direc 'Digital Features 1 UTB\obj_vis_ir_' aux '.mat'],'mask_final');
    %[X_dig,Y_dig]=meshgrid([min(x_dig):0.05:max(x_dig)], [min(y_dig):0.05:max(y_dig)]);
    mask_dig=mask_dig.mask_final;
    
    
    
    
    
    
    X_dig=ones([1024,1280])*NaN;
    Y_dig=ones([1024,1280])*NaN;
    z_est_dig=ones([1024,1280])*NaN;
    t_est_dig=ones([1024,1280])*NaN;
    
    X_dig(mask_dig)=x_dig;
    Y_dig(mask_dig)=y_dig;
    z_est_dig(mask_dig)=z_dig;
    t_est_dig(mask_dig)=t_dig;
    % 
    % [X_dig,Y_dig]=meshgrid(linspace(min(x_dig),max(x_dig),1280), linspace(min(y_dig),max(y_dig),1024));
    % 
    % z_est_dig=griddata(x_dig,y_dig,z_dig,X_dig,Y_dig,'cubic');
    % t_est_dig=griddata(x_dig,y_dig,t_dig,X_dig,Y_dig,'cubic');
    % z_est_dig(z_est_dig>max(z_dig))=NaN;
    % t_est_dig(t_est_dig>max(t_dig))=NaN;
    % z_est_dig(z_est_dig<min(z_dig))=NaN;
    % t_est_dig(t_est_dig<min(t_dig))=NaN;
    
    % Visualización de imagenes 3D y térmicas
    figure;
    subplot(121)
    imagesc(z_est_dig)
    colorbar
    axis equal
    subplot(122)
    imagesc(t_est_dig)
    c = colorbar; 
    c.Label.String = '° C';
    axis equal
    
    % Visualizacion 3D
    figure;
    s = surf(-X_dig,Y_dig,-z_est_dig,'FaceColor', 'interp',...
                'EdgeColor', 'none',...
                'FaceLighting', 'phong');
    
    
    set(gca, 'DataAspectRatio', [1, 1, 1])
    axis equal;
    view(180,90);
    camlight right
    axis on
    grid on
    xlabel('x (mm)')
    ylabel('y (mm)')
    zlabel('z (mm)')
    
    
    
    
    disp('Cantidad de datos NaN:')
    sum(sum(isnan(z_est_dig))) % Ojo con esto
    sum(sum(isnan(t_est_dig))) % Ojo con esto
    
    x_trad=pc_trad.Location(:,1);
    y_trad=pc_trad.Location(:,2);
    z_trad=pc_trad.Location(:,3);
    t_trad=pc_trad.Intensity;
    mask_trad=load([direc 'Two Stages 2\obj_vis_ir_' aux '.mat'],'mask_final');
    %[X_trad,Y_trad]=meshgrid([min(x_trad):0.05:max(x_trad)], [min(y_trad):0.05:max(y_trad)]);
    mask_trad=mask_trad.mask_final;
    
    X_trad=ones([1024,1280])*NaN;
    Y_trad=ones([1024,1280])*NaN;
    z_est_trad=ones([1024,1280])*NaN;
    t_est_trad=ones([1024,1280])*NaN;
    
    X_trad(mask_trad)=x_trad;
    Y_trad(mask_trad)=y_trad;
    z_est_trad(mask_trad)=z_trad;
    t_est_trad(mask_trad)=t_trad;
    % 
    % [X_trad,Y_trad]=meshgrid(linspace(min(x_trad),max(x_trad),1280), linspace(min(y_trad),max(y_trad),1024));
    % 
    % z_est_trad=griddata(x_trad,y_trad,z_trad,X_trad,Y_trad,'cubic');
    % t_est_trad=griddata(x_trad,y_trad,t_trad,X_trad,Y_trad,'cubic');
    % z_est_trad(z_est_trad>max(z_trad))=NaN;
    % t_est_trad(t_est_trad>max(t_trad))=NaN;
    % z_est_trad(z_est_trad<min(z_trad))=NaN;
    % t_est_trad(t_est_trad<min(t_trad))=NaN;
    
    % Visualización de imagenes 3D y térmicas
    figure;
    subplot(121)
    imagesc(z_est_trad)
    colorbar
    axis equal
    subplot(122)
    imagesc(t_est_trad)
    c = colorbar; 
    c.Label.String = '° C';
    axis equal
    
    % Visualizacion 3D
    figure;
    s = surf(-X_trad,Y_trad,-z_est_trad,'FaceColor', 'interp',...
                'EdgeColor', 'none',...
                'FaceLighting', 'phong');
    
    set(gca, 'DataAspectRatio', [1, 1, 1])
    axis equal;
    view(180,90);
    camlight right
    axis on
    grid on
    xlabel('x (mm)')
    ylabel('y (mm)')
    zlabel('z (mm)')
    
    disp('Cantidad de datos NaN:')
    sum(sum(isnan(z_est_trad))) % Ojo con esto
    sum(sum(isnan(t_est_trad))) % Ojo con esto

    sgtitle(['Pose ' aux])


    X_dig(~BW)=NaN;
    Y_dig(~BW)=NaN;
    z_est_dig(~BW)=NaN;
    t_est_dig(~BW)=NaN;




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

    subplot(231)
    s = surf(-X_dig,Y_dig,-z_est_dig,'FaceColor', 'interp',...
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

%     subplot(235)
%     s = surf(map_error,'FaceColor', 'interp',...
%                     'EdgeColor', 'none',...
%                     'FaceLighting', 'phong');
% 
    subplot(233)
    imagesc(map_error)
    c = colorbar; 
    c.Label.String = 'Error (mm)';
    title('Z error')

    %subplot(232)
%     figure(8);
%     %m_mask=reshape(m_mask,[length(tamano),length(tamano)]);
%         m_mask=reshape(m_mask,[1024,1280]);
%     %imshow(m_mask)
%     histogram(error,'EdgeColor','none')
%     disp(['Error dig has ' num2str(size(error))])
%     hold on
%     title('Hist')
    %imhist(map_error)
    sgtitle(['Pose ' num2str(i-1)])
    error_rms(i)=sqrt((sum(error.^2))/length(error));
    disp(['RMSE: ' num2str(error_rms(i))])

    %pause(2)
    %close all


    figure(7);
    % Trad

    X_trad(~BW)=NaN;
    Y_trad(~BW)=NaN;
    z_est_trad(~BW)=NaN;
    t_est_trad(~BW)=NaN;




%    Error
    %tamano=300:700;
    %X_dig=X_dig(tamano,tamano);
    %Y_dig=Y_dig(tamano,tamano);
    %z_est_dig=z_est_dig(tamano,tamano);
    x=X_trad(:);
    y=Y_trad(:);
    z=z_est_trad(:);
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
    error_trad=z_hat-z_n;
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
            error_2d(j)=error_trad(cnt_data);
            m_mask(j)=1;
        end
    end
    %map_error=reshape(error_2d,[length(tamano),length(tamano)]);
    %z_est=reshape(z_est,[length(tamano),length(tamano)]);
        map_error=reshape(error_2d,[1024,1280]);
    z_est=reshape(z_est,[1024,1280]);

    subplot(234)
    s = surf(-X_dig,Y_dig,-z_est_dig,'FaceColor', 'interp',...
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

%     subplot(235)
%     s = surf(map_error,'FaceColor', 'interp',...
%                     'EdgeColor', 'none',...
%                     'FaceLighting', 'phong');
% 
    subplot(236)
    imagesc(map_error)
    c = colorbar; 
    c.Label.String = 'Error (mm)';
    title('Z error')

    %subplot(232)
    figure(8);
    histogram(error,'EdgeColor','none')
    hold on


    %m_mask=reshape(m_mask,[length(tamano),length(tamano)]);
        m_mask=reshape(m_mask,[1024,1280]);
    %imshow(m_mask)
    histogram(error_trad,'EdgeColor','none')
    size(error_trad)
    disp(['Error trad has ' num2str(size(error_trad))])
    hold on
    title('Hist')
    %imhist(map_error)
    sgtitle(['Pose ' num2str(i-1)])
    error_rms(i+1)=sqrt((sum(error_trad.^2))/length(error_trad));
    disp(['RMSE: ' num2str(error_rms(i+1))])
%    pause();
%    close all
%end

% figure;
% plot(error_rms([8:20 22:end]))
% xlabel('Pose')
% ylabel('RMSE')
% title('RMSE')
% 
% mean_error=mean(error_rms([8:20 22:end]));
% desv=std(error_rms([8:20 22:end]));
% disp(['RMSE= ' num2str(mean_error) '±' num2str(desv)])
% 
% figure; 
%     s = surf(map_error,'FaceColor', 'interp',...
%                     'EdgeColor', 'none',...
%                     'FaceLighting', 'phong');


figure(9); clf; hold on
commonEdges = -0.25:0.0025:0.25;           % mismo rango y paso para ambos
histogram(error,      'BinEdges',commonEdges,...
                      'FaceAlpha',0.6,'EdgeColor','none');
histogram(error_trad, 'BinEdges',commonEdges,...
                      'FaceAlpha',0.6,'EdgeColor','none');
legend('Proposed','Conventional')
xlabel('Error (mm)'); ylabel('Counts')
%title(sprintf('Pose %02d – Histograma comparativo',i))
hold off
%%

error_rms=[];
error_rms_trad=[];
for i=8:28

    if i>9
        aux=num2str(i);
    else
        aux=['0' num2str(i)];
    end
    
% Cargamos digital features image
    pc_dig = pcread([direc 'Digital Features 1 UTB\new_obj_vis_ir_' aux '.ply']);
    figure; 
    vis_texture=imread([direc 'pose_' aux '\im_51.png']);
    imshow(vis_texture)
    level = graythresh(vis_texture);        %Umbralizamos
    BW=imbinarize(vis_texture,level);
%     figure;
%     imshow(BW)
    mirror_corners = find_corners_mirror(BW,1280);
    mirror_corners(1,1)=mirror_corners(1,1)+10;
    mirror_corners(1,4)=mirror_corners(1,4)+10;
    mirror_corners(1,2)=mirror_corners(1,2)-10;
    mirror_corners(1,3)=mirror_corners(1,3)-10;
    mirror_corners(2,1)=mirror_corners(2,1)+10;
    mirror_corners(2,4)=mirror_corners(2,4)-10;
    mirror_corners(2,2)=mirror_corners(2,2)+10;
    mirror_corners(2,3)=mirror_corners(2,3)-10;
    hold on
    %plot(mirror_corners(1,:),mirror_corners(2,:),'+r')
    line([mirror_corners(1,:) mirror_corners(1,1)],[mirror_corners(2,:) mirror_corners(2,1)],'Color','blue','LineStyle','--','Linewidth',1.5)
    %BW = zeros([1024,1280]);
    mirror_corners=double(mirror_corners);
    
    % Usar la función poly2mask para crear la máscara a partir de las coordenadas
    BW = poly2mask(mirror_corners(1,:)', mirror_corners(2,:)', 1024, 1280);
    figure;
    imshow(BW)
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
    % Creamos matrices de información
    % Digital features I
    x_dig=pc_dig.Location(:,1);
    y_dig=pc_dig.Location(:,2);
    z_dig=pc_dig.Location(:,3);
    t_dig=pc_dig.Intensity;
    mask_dig=load([direc 'Digital Features 1 UTB\obj_vis_ir_' aux '.mat'],'mask_final');
    %[X_dig,Y_dig]=meshgrid([min(x_dig):0.05:max(x_dig)], [min(y_dig):0.05:max(y_dig)]);
    mask_dig=mask_dig.mask_final;
    
    
    
    
    
    
    X_dig=ones([1024,1280])*NaN;
    Y_dig=ones([1024,1280])*NaN;
    z_est_dig=ones([1024,1280])*NaN;
    t_est_dig=ones([1024,1280])*NaN;
    
    X_dig(mask_dig)=x_dig;
    Y_dig(mask_dig)=y_dig;
    z_est_dig(mask_dig)=z_dig;
    t_est_dig(mask_dig)=t_dig;
    % 
    % [X_dig,Y_dig]=meshgrid(linspace(min(x_dig),max(x_dig),1280), linspace(min(y_dig),max(y_dig),1024));
    % 
    % z_est_dig=griddata(x_dig,y_dig,z_dig,X_dig,Y_dig,'cubic');
    % t_est_dig=griddata(x_dig,y_dig,t_dig,X_dig,Y_dig,'cubic');
    % z_est_dig(z_est_dig>max(z_dig))=NaN;
    % t_est_dig(t_est_dig>max(t_dig))=NaN;
    % z_est_dig(z_est_dig<min(z_dig))=NaN;
    % t_est_dig(t_est_dig<min(t_dig))=NaN;
    
    % Visualización de imagenes 3D y térmicas
    figure;
    subplot(121)
    imagesc(z_est_dig)
    colorbar
    axis equal
    subplot(122)
    imagesc(t_est_dig)
    c = colorbar; 
    c.Label.String = '° C';
    axis equal
    
    % Visualizacion 3D
    figure;
    s = surf(-X_dig,Y_dig,-z_est_dig,'FaceColor', 'interp',...
                'EdgeColor', 'none',...
                'FaceLighting', 'phong');
    
    
    set(gca, 'DataAspectRatio', [1, 1, 1])
    axis equal;
    view(180,90);
    camlight right
    axis on
    grid on
    xlabel('x (mm)')
    ylabel('y (mm)')
    zlabel('z (mm)')
    
    
    
    
    disp('Cantidad de datos NaN:')
    sum(sum(isnan(z_est_dig))) % Ojo con esto
    sum(sum(isnan(t_est_dig))) % Ojo con esto
    
    x_trad=pc_trad.Location(:,1);
    y_trad=pc_trad.Location(:,2);
    z_trad=pc_trad.Location(:,3);
    t_trad=pc_trad.Intensity;
    mask_trad=load([direc 'Two Stages 2\obj_vis_ir_' aux '.mat'],'mask_final');
    %[X_trad,Y_trad]=meshgrid([min(x_trad):0.05:max(x_trad)], [min(y_trad):0.05:max(y_trad)]);
    mask_trad=mask_trad.mask_final;
    
    X_trad=ones([1024,1280])*NaN;
    Y_trad=ones([1024,1280])*NaN;
    z_est_trad=ones([1024,1280])*NaN;
    t_est_trad=ones([1024,1280])*NaN;
    
    X_trad(mask_trad)=x_trad;
    Y_trad(mask_trad)=y_trad;
    z_est_trad(mask_trad)=z_trad;
    t_est_trad(mask_trad)=t_trad;
    % 
    % [X_trad,Y_trad]=meshgrid(linspace(min(x_trad),max(x_trad),1280), linspace(min(y_trad),max(y_trad),1024));
    % 
    % z_est_trad=griddata(x_trad,y_trad,z_trad,X_trad,Y_trad,'cubic');
    % t_est_trad=griddata(x_trad,y_trad,t_trad,X_trad,Y_trad,'cubic');
    % z_est_trad(z_est_trad>max(z_trad))=NaN;
    % t_est_trad(t_est_trad>max(t_trad))=NaN;
    % z_est_trad(z_est_trad<min(z_trad))=NaN;
    % t_est_trad(t_est_trad<min(t_trad))=NaN;
    
    % Visualización de imagenes 3D y térmicas
    figure;
    subplot(121)
    imagesc(z_est_trad)
    colorbar
    axis equal
    subplot(122)
    imagesc(t_est_trad)
    c = colorbar; 
    c.Label.String = '° C';
    axis equal
    
    % Visualizacion 3D
    figure;
    s = surf(-X_trad,Y_trad,-z_est_trad,'FaceColor', 'interp',...
                'EdgeColor', 'none',...
                'FaceLighting', 'phong');
    
    set(gca, 'DataAspectRatio', [1, 1, 1])
    axis equal;
    view(180,90);
    camlight right
    axis on
    grid on
    xlabel('x (mm)')
    ylabel('y (mm)')
    zlabel('z (mm)')
    
    disp('Cantidad de datos NaN:')
    sum(sum(isnan(z_est_trad))) % Ojo con esto
    sum(sum(isnan(t_est_trad))) % Ojo con esto

    sgtitle(['Pose ' aux])


    X_dig(~BW)=NaN;
    Y_dig(~BW)=NaN;
    z_est_dig(~BW)=NaN;
    t_est_dig(~BW)=NaN;




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

    subplot(231)
    s = surf(-X_dig,Y_dig,-z_est_dig,'FaceColor', 'interp',...
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

%     subplot(235)
%     s = surf(map_error,'FaceColor', 'interp',...
%                     'EdgeColor', 'none',...
%                     'FaceLighting', 'phong');
% 
    subplot(233)
    imagesc(map_error)
    c = colorbar; 
    c.Label.String = 'Error (mm)';
    title('Z error')

    %subplot(232)
%     figure(8);
%     %m_mask=reshape(m_mask,[length(tamano),length(tamano)]);
%         m_mask=reshape(m_mask,[1024,1280]);
%     %imshow(m_mask)
%     histogram(error,'EdgeColor','none')
%     disp(['Error dig has ' num2str(size(error))])
%     hold on
%     title('Hist')
    %imhist(map_error)
    sgtitle(['Pose ' num2str(i-1)])
    error_rms(i)=sqrt((sum(error.^2))/length(error));
    disp(['RMSE: ' num2str(error_rms(i))])

    %pause(2)
    %close all



    % Trad

    X_trad(~BW)=NaN;
    Y_trad(~BW)=NaN;
    z_est_trad(~BW)=NaN;
    t_est_trad(~BW)=NaN;




%    Error
    %tamano=300:700;
    %X_dig=X_dig(tamano,tamano);
    %Y_dig=Y_dig(tamano,tamano);
    %z_est_dig=z_est_dig(tamano,tamano);
    x=X_trad(:);
    y=Y_trad(:);
    z=z_est_trad(:);
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
    error_trad=z_hat-z_n;
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
            error_2d(j)=error_trad(cnt_data);
            m_mask(j)=1;
        end
    end
    %map_error=reshape(error_2d,[length(tamano),length(tamano)]);
    %z_est=reshape(z_est,[length(tamano),length(tamano)]);
        map_error=reshape(error_2d,[1024,1280]);
    z_est=reshape(z_est,[1024,1280]);

    subplot(234)
    s = surf(-X_dig,Y_dig,-z_est_dig,'FaceColor', 'interp',...
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

%     subplot(235)
%     s = surf(map_error,'FaceColor', 'interp',...
%                     'EdgeColor', 'none',...
%                     'FaceLighting', 'phong');
% 
    subplot(236)
    imagesc(map_error)
    c = colorbar; 
    c.Label.String = 'Error (mm)';
    title('Z error')

    %subplot(232)
%     figure(8);
%     histogram(error,'EdgeColor','none')
%     hold on


    %m_mask=reshape(m_mask,[length(tamano),length(tamano)]);
        m_mask=reshape(m_mask,[1024,1280]);
    %imshow(m_mask)
%     histogram(error_trad,'EdgeColor','none')
%     size(error_trad)
     disp(['Error trad has ' num2str(size(error_trad))])
%     hold on
%     title('Hist')
    %imhist(map_error)
    sgtitle(['Pose ' num2str(i-1)])
    error_rms_trad(i)=sqrt((sum(error_trad.^2))/length(error_trad));
    disp(['RMSE: ' num2str(error_rms_trad(i))])
%    pause();
%    close all
%end

% figure;
% plot(error_rms([8:20 22:end]))
% xlabel('Pose')
% ylabel('RMSE')
% title('RMSE')
% 
% mean_error=mean(error_rms([8:20 22:end]));
% desv=std(error_rms([8:20 22:end]));
% disp(['RMSE= ' num2str(mean_error) '±' num2str(desv)])
% 
% figure; 
%     s = surf(map_error,'FaceColor', 'interp',...
%                     'EdgeColor', 'none',...
%                     'FaceLighting', 'phong');


%figure(9); clf; hold on
    subplot(235)
    commonEdges = -0.20:0.0025:0.20;           % mismo rango y paso para ambos
    histogram(error,      'BinEdges',commonEdges,...
                          'FaceAlpha',0.6,'EdgeColor','none');
    hold on
    histogram(error_trad, 'BinEdges',commonEdges,...
                          'FaceAlpha',0.6,'EdgeColor','none');
    legend('Proposed','Conventional')
    xlabel('Error (mm)'); ylabel('Counts')
    %title(sprintf('Pose %02d – Histograma comparativo',i))
    hold off
    pause(0.1);
    close all

end
error_rms_trad(21)=[];
error_rms(21)=[];
error_rms_trad(1:7)=[];
error_rms(1:7)=[];
figure;
plot(error_rms,'b','LineWidth',1.5)
hold on
plot(error_rms_trad,'r','LineWidth',1.5)
ylabel('RMS error (mm)')
xlabel('Poses')
legend({'Proposed','Conventional'})