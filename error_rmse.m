close all
clc
clear all
direc='D:\R3D\Adquisiciones\Experimento 54\Map_cal\SL\';
n_poses=62;
addpath 'C:\Users\Eberto Benjumea\Desktop\OneDrive - Universidad Tecnológica de Bolívar\Doctorado\Proyectos\Codigos Reconstruccion 3D Learning\Codigos Tesis All Digital Features'

% %%
% figure;
% for i=1:n_poses
%     if i-1<10
%         aux=['0' num2str(i-1)];
%     else
%         aux=num2str(i-1);
%     end
%     disp(['Pose ' aux])
%     %subplot(121)
%     % Cargamos digital features image
%     pc_dig = pcread([direc 'Digital Features 1 UTB\obj_vis_ir_' aux '.ply']);
% %     pcshow(pc_dig);
% %     xlabel('x (mm)')
% %     ylabel('y (mm)')
% %     zlabel('z (mm)')
% %     c = colorbar; 
% %     c.Label.String = '° C';
% %     set(gca, 'Zdir', 'reverse')
% %     set(gca, 'Xdir', 'reverse')
% %     view(30,45)
% %     grid off
% %     zlim([350 430])
% %     title('Digital features I')
%     
%     % Cargamos traditional image
%     %subplot(122)
%     pc_trad = pcread([direc 'Two Stages 2\obj_vis_ir_' aux '.ply']);
% %     pcshow(pc_trad);
% %     xlabel('x (mm)')
% %     ylabel('y (mm)')
% %     zlabel('z (mm)')
% %     c = colorbar; 
% %     c.Label.String = '° C';
% %     set(gca, 'Zdir', 'reverse')
% %     set(gca, 'Xdir', 'reverse')
% %     view(30,45)
% %     grid off
% %     zlim([350 430])
% %     title('Traditional')
%     
%     % Digital features I
%     x_dig=pc_dig.Location(:,1);
%     y_dig=pc_dig.Location(:,2);
%     z_dig=pc_dig.Location(:,3);
%     t_dig=pc_dig.Intensity;
%     %[X_dig,Y_dig]=meshgrid([min(x_dig):0.05:max(x_dig)], [min(y_dig):0.05:max(y_dig)]);
%     
%     [X_dig,Y_dig]=meshgrid(linspace(min(x_dig),max(x_dig),1280), linspace(min(y_dig),max(y_dig),1024));
%     
%     
%     
%     z_est_dig=griddata(x_dig,y_dig,z_dig,X_dig,Y_dig,'cubic');
%     t_est_dig=griddata(x_dig,y_dig,t_dig,X_dig,Y_dig,'cubic');
%     z_est_dig(z_est_dig>max(z_dig))=NaN;
%     t_est_dig(t_est_dig>max(t_dig))=NaN;
%     z_est_dig(z_est_dig<min(z_dig))=NaN;
%     t_est_dig(t_est_dig<min(t_dig))=NaN;
%     
% 
%     
%     % Visualizacion 3D
%     subplot(121)
%     s = surf(-X_dig,Y_dig,-z_est_dig,'FaceColor', 'interp',...
%                 'EdgeColor', 'none',...
%                 'FaceLighting', 'phong');
%     
%     
%     set(gca, 'DataAspectRatio', [1, 1, 1])
%     axis equal;
%     view(180,90);
%     camlight right
%     axis on
%     grid on
%     xlabel('x (mm)')
%     ylabel('y (mm)')
%     zlabel('z (mm)')
%     
%     % Traditional
%     x_trad=pc_trad.Location(:,1);
%     y_trad=pc_trad.Location(:,2);
%     z_trad=pc_trad.Location(:,3);
%     t_trad=pc_trad.Intensity;
%     %[X_trad,Y_trad]=meshgrid([min(x_trad):0.05:max(x_trad)], [min(y_trad):0.05:max(y_trad)]);
%     [X_trad,Y_trad]=meshgrid(linspace(min(x_trad),max(x_trad),1280), linspace(min(y_trad),max(y_trad),1024));
%     
%     
%     
%     z_est_trad=griddata(x_trad,y_trad,z_trad,X_trad,Y_trad,'cubic');
%     t_est_trad=griddata(x_trad,y_trad,t_trad,X_trad,Y_trad,'cubic');
%     z_est_trad(z_est_trad>max(z_trad))=NaN;
%     t_est_trad(t_est_trad>max(t_trad))=NaN;
%     z_est_trad(z_est_trad<min(z_trad))=NaN;
%     t_est_trad(t_est_trad<min(t_trad))=NaN;
%     
%     % Visualizacion 3D
%     subplot(122)
%     s = surf(-X_trad,Y_trad,-z_est_trad,'FaceColor', 'interp',...
%                 'EdgeColor', 'none',...
%                 'FaceLighting', 'phong');
%     
%     
%     set(gca, 'DataAspectRatio', [1, 1, 1])
%     axis equal;
%     view(180,90);
%     camlight right
%     axis on
%     grid on
%     xlabel('x (mm)')
%     ylabel('y (mm)')
%     zlabel('z (mm)')
% 
% 
% 
%     pause();
% end
% 
% %%
% close all
% poses=[5 8 13 14 16 17 18 19 20 21 24 25 26 28 29 30 31 32 33 35 36 39 40 41 42 43 44 46 47 49 53 60];
% length(poses)
% figure;
% for i=1:length(poses)
%     if poses(i)<10
%         aux=['0' num2str(poses(i))];
%     else
%         aux=num2str(poses(i));
%     end
%     disp(['Pose ' aux])
%     
%     subplot(121)
%     % Cargamos digital features image
%     pc_dig = pcread([direc 'Digital Features 1 UTB\new_obj_vis_ir_' aux '.ply']);
%     pcshow(pc_dig);
%     xlabel('x (mm)')
%     ylabel('y (mm)')
%     zlabel('z (mm)')
%     c = colorbar; 
%     c.Label.String = '° C';
%     set(gca, 'Zdir', 'reverse')
%     set(gca, 'Xdir', 'reverse')
%     view(30,45)
%     grid off
%     zlim([350 430])
%     title('Digital features I')
%     
%     % Cargamos traditional image
%     subplot(122)
%     pc_trad = pcread([direc 'Two Stages 2\new_obj_vis_ir_' aux '.ply']);
%     pcshow(pc_trad);
%     xlabel('x (mm)')
%     ylabel('y (mm)')
%     zlabel('z (mm)')
%     c = colorbar; 
%     c.Label.String = '° C';
%     set(gca, 'Zdir', 'reverse')
%     set(gca, 'Xdir', 'reverse')
%     view(30,45)
%     grid off
%     zlim([350 430])
%     title('Traditional')
%     set(gcf,'Color','white')
%     sgtitle(['Pose ' aux])
% %     % Digital features I
% %     x_dig=pc_dig.Location(:,1);
% %     y_dig=pc_dig.Location(:,2);
% %     z_dig=pc_dig.Location(:,3);
% %     t_dig=pc_dig.Intensity;
% %     %[X_dig,Y_dig]=meshgrid([min(x_dig):0.05:max(x_dig)], [min(y_dig):0.05:max(y_dig)]);
% %     
% %     [X_dig,Y_dig]=meshgrid(linspace(min(x_dig),max(x_dig),1280), linspace(min(y_dig),max(y_dig),1024));
% %     
% %     
% %     
% %     z_est_dig=griddata(x_dig,y_dig,z_dig,X_dig,Y_dig,'cubic');
% %     t_est_dig=griddata(x_dig,y_dig,t_dig,X_dig,Y_dig,'cubic');
% %     z_est_dig(z_est_dig>max(z_dig))=NaN;
% %     t_est_dig(t_est_dig>max(t_dig))=NaN;
% %     z_est_dig(z_est_dig<min(z_dig))=NaN;
% %     t_est_dig(t_est_dig<min(t_dig))=NaN;
% %     
% % 
% %     
% %     % Visualizacion 3D
% %     subplot(121)
% %     s = surf(-X_dig,Y_dig,-z_est_dig,'FaceColor', 'interp',...
% %                 'EdgeColor', 'none',...
% %                 'FaceLighting', 'phong');
% %     
% %     
% %     set(gca, 'DataAspectRatio', [1, 1, 1])
% %     axis equal;
% %     view(180,90);
% %     camlight right
% %     axis on
% %     grid on
% %     xlabel('x (mm)')
% %     ylabel('y (mm)')
% %     zlabel('z (mm)')
% %     
% %     % Traditional
% %     x_trad=pc_trad.Location(:,1);
% %     y_trad=pc_trad.Location(:,2);
% %     z_trad=pc_trad.Location(:,3);
% %     t_trad=pc_trad.Intensity;
% %     %[X_trad,Y_trad]=meshgrid([min(x_trad):0.05:max(x_trad)], [min(y_trad):0.05:max(y_trad)]);
% %     [X_trad,Y_trad]=meshgrid(linspace(min(x_trad),max(x_trad),1280), linspace(min(y_trad),max(y_trad),1024));
% %     
% %     
% %     
% %     z_est_trad=griddata(x_trad,y_trad,z_trad,X_trad,Y_trad,'cubic');
% %     t_est_trad=griddata(x_trad,y_trad,t_trad,X_trad,Y_trad,'cubic');
% %     z_est_trad(z_est_trad>max(z_trad))=NaN;
% %     t_est_trad(t_est_trad>max(t_trad))=NaN;
% %     z_est_trad(z_est_trad<min(z_trad))=NaN;
% %     t_est_trad(t_est_trad<min(t_trad))=NaN;
% %     
% %     % Visualizacion 3D
% %     subplot(122)
% %     s = surf(-X_trad,Y_trad,-z_est_trad,'FaceColor', 'interp',...
% %                 'EdgeColor', 'none',...
% %                 'FaceLighting', 'phong');
% %     
% %     
% %     set(gca, 'DataAspectRatio', [1, 1, 1])
% %     axis equal;
% %     view(180,90);
% %     camlight right
% %     axis on
% %     grid on
% %     xlabel('x (mm)')
% %     ylabel('y (mm)')
% %     zlabel('z (mm)')
% 
% 
% 
%     pause(1);
% end

%%
close all
poses=[5 8 18 19 20 24 28 29 30 33 41 42 47];
poses=[5 8 18 19 20 24 28 29 30 41 42 47];
length(poses)
X_digcell={};
Y_digcell={};
z_est_digcell={};
t_est_digcell={};
X_tradcell={};
Y_tradcell={};
z_est_tradcell={};
t_est_tradcell={};
for i=1:length(poses)
    if poses(i)<10
        aux=['0' num2str(poses(i))];
    else
        aux=num2str(poses(i));
    end
    disp(['Pose ' aux])
    

    % Cargamos digital features image
    pc_dig = pcread([direc 'Digital Features 1 UTB\new_obj_vis_ir_' aux '.ply']);
    % Cargamos traditional image
    pc_trad = pcread([direc 'Two Stages 2\new_obj_vis_ir_' aux '.ply']);

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
    

    
    % Visualizacion 3D
    subplot(221)
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
    title('3D data by digital features')

    subplot(222)
    imagesc(t_est_dig)
    c = colorbar;
    c.Label.String = '° C';
    title('t by digital features')
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
    subplot(223)
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
    title('3D data by traditional method')

    subplot(224)
    imagesc(t_est_trad)
    c = colorbar;
    c.Label.String = '° C';
    title('t by traditional method')

    sgtitle(['Pose ' aux])

    X_digcell{i}=X_dig;
    Y_digcell{i}=Y_dig;
    z_est_digcell{i}=z_est_dig;
    t_est_digcell{i}=t_est_dig;
    X_tradcell{i}=X_trad;
    Y_tradcell{i}=Y_trad;
    z_est_tradcell{i}=z_est_trad;
    t_est_tradcell{i}=t_est_trad;

    pause(5);
end


for i=1:length(z_est_digcell)

    mean_dig=mean(mean(z_est_digcell{i},'omitnan'),'omitnan');
    mean_trad=mean(mean(z_est_tradcell{i},'omitnan'),'omitnan');

    z_dig_norm=z_est_digcell{i}-mean_dig;
    z_trad_norm=z_est_tradcell{i}-mean_trad;

    sum(sum(isnan(z_dig_norm)))
    sum(sum(isnan(z_trad_norm)))

    %figure, imshow(mask_dig)
    difference=z_dig_norm-z_trad_norm;
    mask_dig=isnan(difference);
    dif=(sum(sum((difference).^2,"omitnan"),"omitnan"))/length(difference(~isnan(difference)));
    rmse_z(i)=sqrt(dif);



    mean_dig=mean(mean(t_est_digcell{i},'omitnan'),'omitnan');
    mean_trad=mean(mean(t_est_tradcell{i},'omitnan'),'omitnan');

    t_dig_norm=t_est_digcell{i}-mean_dig;
    t_trad_norm=t_est_tradcell{i}-mean_trad;

    sum(sum(isnan(t_dig_norm)))
    sum(sum(isnan(t_trad_norm)))

    %figure, imshow(mask_dig)
    difference=t_dig_norm-t_trad_norm;
    mask_dig=isnan(difference);
    dif=(sum(sum((difference).^2,"omitnan"),"omitnan"))/length(difference(~isnan(difference)));
    rmse_t(i)=sqrt(dif);

%     figure;
%     %subplot(121)
%         s = surf(-X_tradcell{i},Y_tradcell{i},-z_trad_norm,'FaceColor', 'b',...
%                 'EdgeColor', 'none',...
%                 'FaceLighting', 'phong');
%     
%     
%     set(gca, 'DataAspectRatio', [1, 1, 1])
%     axis equal;
%     view(180,90);
%     camlight right
%     axis on
%     grid on
%     xlabel('x (mm)')
%     ylabel('y (mm)')
%     zlabel('z (mm)')
%     title('3D data by traditional method')
%     hold on
%     %subplot(122)
%     s = surf(-X_digcell{i},Y_digcell{i},-z_dig_norm,'FaceColor', 'r',...
%                 'EdgeColor', 'none',...
%                 'FaceLighting', 'phong');
%     
%     
%     set(gca, 'DataAspectRatio', [1, 1, 1])
%     axis equal;
%     view(180,90);
%     camlight right
%     axis on
%     grid on
%     xlabel('x (mm)')
%     ylabel('y (mm)')
%     zlabel('z (mm)')
%     title('3D data by digital features')


end

figure;
for i=1:length(z_est_digcell)
    if poses(i)<10
        aux=['0' num2str(poses(i))];
    else
        aux=num2str(poses(i));
    end
%     subplot(121)
%     pc_dig = pcread([direc 'Digital Features 1 UTB\new_obj_vis_ir_' aux '.ply']);
%     pcshow(pc_dig);
%     xlabel('x (mm)')
%     ylabel('y (mm)')
%     zlabel('z (mm)')
%     c = colorbar; 
%     c.Label.String = '° C';
%     set(gca, 'Zdir', 'reverse')
%     set(gca, 'Xdir', 'reverse')
%     view(30,45)
%     grid off
%     zlim([350 430])
%     title('Digital')
%     set(gcf,'Color','white')
% 
%     subplot(122)
%     plot(rmse_z(1:i),'-or')
%     xlim([0 16])

    pause(1);






    subplot(1, 2, 1) % top row of the 3x3 grid
    pc_dig = pcread([direc 'Digital Features 1 UTB\new_obj_vis_ir_' aux '.ply']);
    pcshow(pc_dig,'BackgroundColor','w');
    xlabel('x (mm)')
    ylabel('y (mm)')
    zlabel('z (mm)')
    c = colorbar; 
    c.Label.String = '° C';
    set(gca, 'Zdir', 'reverse')
    set(gca, 'Xdir', 'reverse')
    view(180,90)
    %view(30,45)
    grid off
    zlim([350 430])
    title('Digital')
    %set(gcf,'Color','white')
%     subplot(1, 2, 2) % bottom 2/3 of the first column of the grid
%     plot(rmse_z(1:i),'-or')
%     xlabel('Poses')
%     ylabel('RMSE (mm)')
%     xlim([0 16])
    subplot(1, 2, 2); % the rest of the grid
    plot(rmse_t(1:i),'-ob')
    xlim([0 14])
    xlabel('Poses')
    ylabel('RMSE (° C)')
end
subplot(1,2,2); 

yl = yline(mean(rmse_t),'--',{'Mean: ', [num2str(round(mean(rmse_t),2)) ' mm.']},'LabelVerticalAlignment','top','LineWidth',1);
yl.FontSize = 7;

%% Solicitado Prof. Marrugo - Error RMS en x
close all
figure;
detrend_z_est_digcell={};
detrend_t_est_digcell={};
for i=1:length(z_est_digcell)
    imagesc(z_est_digcell{i})
    for j=1:size(z_est_digcell{i},1)
        mean_dig=mean(z_est_digcell{i}(j,:),'omitnan');
        aux_dig=z_est_digcell{i}(j,:)-mean_dig;
        [z_est_final,t_est_final]=detrend_profile_x(aux_dig, t_est_digcell{i}(j,:));
        detrend_z_est_digcell{i}(j,:)=z_est_final;
        detrend_t_est_digcell{i}(j,:)=t_est_final;        

    end
    imagesc(detrend_z_est_digcell{i})
    title(['Pose ' num2str(i)])
    pause(0.5);
end

figure;
detrend_z_est_tradcell={};
detrend_t_est_tradcell={};
for i=1:length(z_est_tradcell)
    imagesc(z_est_tradcell{i})
    for j=1:size(z_est_tradcell{i},1)
        mean_trad=mean(z_est_tradcell{i}(j,:),'omitnan');
        aux_trad=z_est_tradcell{i}(j,:)-mean_trad;
        [z_est_final,t_est_final]=detrend_profile_x(aux_trad, t_est_tradcell{i}(j,:));
        detrend_z_est_tradcell{i}(j,:)=z_est_final;
        detrend_t_est_tradcell{i}(j,:)=t_est_final;        

    end
    imagesc(detrend_z_est_tradcell{i})
    title(['Pose ' num2str(i)])
    pause(0.5);
end

for i=1:length(detrend_z_est_tradcell)

    mean_dig=mean(mean(detrend_z_est_digcell{i},'omitnan'),'omitnan');
    mean_trad=mean(mean(detrend_z_est_tradcell{i},'omitnan'),'omitnan');

    z_dig_norm=detrend_z_est_digcell{i};%-mean_dig;
    z_trad_norm=detrend_z_est_tradcell{i};%-mean_trad;



    %figure, imshow(mask_dig)
    difference=z_dig_norm-z_trad_norm;
    mask_dig=isnan(difference);
    dif=(sum(sum((difference).^2,"omitnan"),"omitnan"))/length(difference(~isnan(difference)));
    rmse_z(i)=sqrt(dif);



    mean_dig=mean(mean(detrend_t_est_digcell{i},'omitnan'),'omitnan');
    mean_trad=mean(mean(detrend_t_est_tradcell{i},'omitnan'),'omitnan');

    t_dig_norm=detrend_t_est_digcell{i};%-mean_dig;
    t_trad_norm=detrend_t_est_tradcell{i};%-mean_trad;

    sum(sum(isnan(t_dig_norm)))
    sum(sum(isnan(t_trad_norm)))

    %figure, imshow(mask_dig)
    difference=t_dig_norm-t_trad_norm;
    mask_dig=isnan(difference);
    dif=(sum(sum((difference).^2,"omitnan"),"omitnan"))/length(difference(~isnan(difference)));
    rmse_t(i)=sqrt(dif);

%     figure;
%     %subplot(121)
%         s = surf(-X_tradcell{i},Y_tradcell{i},-z_trad_norm,'FaceColor', 'b',...
%                 'EdgeColor', 'none',...
%                 'FaceLighting', 'phong');
%     
%     
%     set(gca, 'DataAspectRatio', [1, 1, 1])
%     axis equal;
%     view(180,90);
%     camlight right
%     axis on
%     grid on
%     xlabel('x (mm)')
%     ylabel('y (mm)')
%     zlabel('z (mm)')
%     title('3D data by traditional method')
%     hold on
%     %subplot(122)
%     s = surf(-X_digcell{i},Y_digcell{i},-z_dig_norm,'FaceColor', 'r',...
%                 'EdgeColor', 'none',...
%                 'FaceLighting', 'phong');
%     
%     
%     set(gca, 'DataAspectRatio', [1, 1, 1])
%     axis equal;
%     view(180,90);
%     camlight right
%     axis on
%     grid on
%     xlabel('x (mm)')
%     ylabel('y (mm)')
%     zlabel('z (mm)')
%     title('3D data by digital features')


end

folder='C:\Users\Eberto Benjumea\Desktop\OneDrive - Universidad Tecnológica de Bolívar\Images for Events and Papers\Applied Optics\new_figures\gif rmse';
cnt=0;
figure();
for i=1:length(detrend_z_est_digcell)
    if poses(i)<10
        aux=['0' num2str(poses(i))];
    else
        aux=num2str(poses(i));
    end
%     subplot(121)
%     pc_dig = pcread([direc 'Digital Features 1 UTB\new_obj_vis_ir_' aux '.ply']);
%     pcshow(pc_dig);
%     xlabel('x (mm)')
%     ylabel('y (mm)')
%     zlabel('z (mm)')
%     c = colorbar; 
%     c.Label.String = '° C';
%     set(gca, 'Zdir', 'reverse')
%     set(gca, 'Xdir', 'reverse')
%     view(30,45)
%     grid off
%     zlim([350 430])
%     title('Digital')
%     set(gcf,'Color','white')
% 
%     subplot(122)
%     plot(rmse_z(1:i),'-or')
%     xlim([0 16])

    pause(1);




    figure(30);
    %subplot(1, 2, 1) % top row of the 3x3 grid
    % subplot(2, 2, [1 3]) % top row of the 3x3 grid
    pc_dig = pcread([direc 'Digital Features 1 UTB\new_obj_vis_ir_' aux '.ply']);
    pcshow(pc_dig,'BackgroundColor','w');
    xlabel('x (mm)')
    ylabel('y (mm)')
    zlabel('z (mm)')
    c = colorbar; 
    c.Label.String = '° C';
    set(gca, 'Zdir', 'reverse')
    set(gca, 'Xdir', 'reverse')
    view(30,45)
    grid off
    zlim([350 430])
    xlim([-140 160])
    ylim([-140 200])
    %title('Digital')
    %set(gcf,'Color','white')
    figure(31);
    %subplot(1, 2, 2)
    %subplot(2, 2, 2) % bottom 2/3 of the first column of the grid
    plot(rmse_z(1:i),'-ob')
    xlabel('Poses')
    ylabel('RMSE (mm)')
    ylim([0 0.18])
    xlim([0 13])
%     subplot(2, 2, 4); % the rest of the grid
%     plot(rmse_t(1:i),'-ob')
%     xlim([0 13])
%     xlabel('Poses')
%     ylabel('RMSE (° C)')
    cnt=cnt+1;

    F = getframe(gcf);
    if cnt<10
        imwrite(F.cdata,[folder '\im0' num2str(cnt) '.png'])
        pause(0.1)
    else
        imwrite(F.cdata,[folder '\im' num2str(cnt) '.png'])
        pause(0.1)
    end
end
figure(31); %subplot(1,2,2); 

yl = yline(mean(rmse_z),'--',{'Mean: ', [num2str(round(mean(rmse_z),2)) ' mm.']},'LabelVerticalAlignment','top','LineWidth',1);
yl.FontSize = 7;
cnt=cnt+1;

F = getframe(gcf);
if cnt<10
    imwrite(F.cdata,[folder '\im0' num2str(cnt) '.png'])
    pause(0.1)
else
    imwrite(F.cdata,[folder '\im' num2str(cnt) '.png'])
    pause(0.1)
end
% subplot(2,2,4); 
% yl = yline(mean(rmse_t),'--',['Mean: ' num2str(round(mean(rmse_t),2)) ' ° C.'],'LabelVerticalAlignment','top','LineWidth',1);
% yl.FontSize = 7;







%%


%%


%% Solicitado Prof. Marrugo - Error RMS en y
close all
figure;
detrend_z_est_digcell={};
detrend_t_est_digcell={};
for i=1:length(z_est_digcell)
    imagesc(z_est_digcell{i})
    for j=1:size(z_est_digcell{i},2)
        mean_dig=mean(z_est_digcell{i}(:,j),'omitnan');
        aux_dig=z_est_digcell{i}(:,j)-mean_dig;
        [z_est_final,t_est_final]=detrend_profile_y(aux_dig, t_est_digcell{i}(:,j));
        detrend_z_est_digcell{i}(:,j)=z_est_final;
        detrend_t_est_digcell{i}(:,j)=t_est_final;        

    end
    imagesc(detrend_z_est_digcell{i})
    title(['Pose ' num2str(i)])
    pause(0.5);
end

figure;
detrend_z_est_tradcell={};
detrend_t_est_tradcell={};
for i=1:length(z_est_tradcell)
    imagesc(z_est_tradcell{i})
    for j=1:size(z_est_tradcell{i},2)
        mean_trad=mean(z_est_tradcell{i}(:,j),'omitnan');
        aux_trad=z_est_tradcell{i}(:,j)-mean_trad;
        [z_est_final,t_est_final]=detrend_profile_y(aux_trad, t_est_tradcell{i}(:,j));
        detrend_z_est_tradcell{i}(:,j)=z_est_final;
        detrend_t_est_tradcell{i}(:,j)=t_est_final;        

    end
    imagesc(detrend_z_est_tradcell{i})
    title(['Pose ' num2str(i)])
    pause(0.5);
end

for i=1:length(detrend_z_est_tradcell)

    mean_dig=mean(mean(detrend_z_est_digcell{i},'omitnan'),'omitnan');
    mean_trad=mean(mean(detrend_z_est_tradcell{i},'omitnan'),'omitnan');

    z_dig_norm=detrend_z_est_digcell{i};%-mean_dig;
    z_trad_norm=detrend_z_est_tradcell{i};%-mean_trad;



    %figure, imshow(mask_dig)
    difference=z_dig_norm-z_trad_norm;
    mask_dig=isnan(difference);
    dif=(sum(sum((difference).^2,"omitnan"),"omitnan"))/length(difference(~isnan(difference)));
    rmse_z(i)=sqrt(dif);



    mean_dig=mean(mean(detrend_t_est_digcell{i},'omitnan'),'omitnan');
    mean_trad=mean(mean(detrend_t_est_tradcell{i},'omitnan'),'omitnan');

    t_dig_norm=detrend_t_est_digcell{i};%-mean_dig;
    t_trad_norm=detrend_t_est_tradcell{i};%-mean_trad;

    sum(sum(isnan(t_dig_norm)))
    sum(sum(isnan(t_trad_norm)))

    %figure, imshow(mask_dig)
    difference=t_dig_norm-t_trad_norm;
    mask_dig=isnan(difference);
    dif=(sum(sum((difference).^2,"omitnan"),"omitnan"))/length(difference(~isnan(difference)));
    rmse_t(i)=sqrt(dif);

%     figure;
%     %subplot(121)
%         s = surf(-X_tradcell{i},Y_tradcell{i},-z_trad_norm,'FaceColor', 'b',...
%                 'EdgeColor', 'none',...
%                 'FaceLighting', 'phong');
%     
%     
%     set(gca, 'DataAspectRatio', [1, 1, 1])
%     axis equal;
%     view(180,90);
%     camlight right
%     axis on
%     grid on
%     xlabel('x (mm)')
%     ylabel('y (mm)')
%     zlabel('z (mm)')
%     title('3D data by traditional method')
%     hold on
%     %subplot(122)
%     s = surf(-X_digcell{i},Y_digcell{i},-z_dig_norm,'FaceColor', 'r',...
%                 'EdgeColor', 'none',...
%                 'FaceLighting', 'phong');
%     
%     
%     set(gca, 'DataAspectRatio', [1, 1, 1])
%     axis equal;
%     view(180,90);
%     camlight right
%     axis on
%     grid on
%     xlabel('x (mm)')
%     ylabel('y (mm)')
%     zlabel('z (mm)')
%     title('3D data by digital features')


end


figure;
for i=1:length(detrend_z_est_digcell)
    if poses(i)<10
        aux=['0' num2str(poses(i))];
    else
        aux=num2str(poses(i));
    end
%     subplot(121)
%     pc_dig = pcread([direc 'Digital Features 1 UTB\new_obj_vis_ir_' aux '.ply']);
%     pcshow(pc_dig);
%     xlabel('x (mm)')
%     ylabel('y (mm)')
%     zlabel('z (mm)')
%     c = colorbar; 
%     c.Label.String = '° C';
%     set(gca, 'Zdir', 'reverse')
%     set(gca, 'Xdir', 'reverse')
%     view(30,45)
%     grid off
%     zlim([350 430])
%     title('Digital')
%     set(gcf,'Color','white')
% 
%     subplot(122)
%     plot(rmse_z(1:i),'-or')
%     xlim([0 16])

    pause(1);






    subplot(2, 2, [1 3]) % top row of the 3x3 grid
    pc_dig = pcread([direc 'Digital Features 1 UTB\new_obj_vis_ir_' aux '.ply']);
    pcshow(pc_dig,'BackgroundColor','w');
    xlabel('x (mm)')
    ylabel('y (mm)')
    zlabel('z (mm)')
    c = colorbar; 
    c.Label.String = '° C';
    set(gca, 'Zdir', 'reverse')
    set(gca, 'Xdir', 'reverse')
    view(30,45)
    grid off
    zlim([350 430])
    title('Digital')
    %set(gcf,'Color','white')
    subplot(2, 2, 2) % bottom 2/3 of the first column of the grid
    plot(rmse_z(1:i),'-or')
    xlabel('Poses')
    ylabel('RMSE (mm)')
    xlim([0 13])
    subplot(2, 2, 4); % the rest of the grid
    plot(rmse_t(1:i),'-ob')
    xlim([0 13])
    xlabel('Poses')
    ylabel('RMSE (° C)')
end
subplot(2,2,2); 

yl = yline(mean(rmse_z),'--',{'Mean: ', [num2str(round(mean(rmse_z),2)) ' mm.']},'LabelVerticalAlignment','top','LineWidth',1);
yl.FontSize = 7;
subplot(2,2,4); 
yl = yline(mean(rmse_t),'--',['Mean: ' num2str(round(mean(rmse_t),2)) ' ° C.'],'LabelVerticalAlignment','top','LineWidth',1);
yl.FontSize = 7;



%% Error Mapeo
close all
figure(80);


rms_t=[];
for i=1:length(detrend_t_est_tradcell)
    % Detectamos en tradicional
    detrend_t_est_tradcell{i};
    t_est_trad_copy=detrend_t_est_tradcell{i};
    % figure; 
    % imshow(detrend_t_est_tradcell{i},[])
    B = sort(detrend_t_est_tradcell{i}(:));
    C = unique(B);
    detrend_t_est_tradcell{i}(isnan(detrend_t_est_tradcell{i}))=C(2); %mean(mean(t_sin_perfil,"omitnan"),"omitnan");
    detrend_t_est_tradcell{i}(detrend_t_est_tradcell{i}==0)=C(2); %mean(mean(t_sin_perfil,"omitnan"),"omitnan");
    
    X=detrend_t_est_tradcell{i};
    Xmin = min(X(:));
    Xmax = max(X(:));
    if isequal(Xmax,Xmin)
        X = 0*X;
    else
        X = (X - Xmin) ./ (Xmax - Xmin);
    end
    
    % Threshold image - manual threshold
    BW = X > 4.588200e-01;
    
    % Open mask with disk
    radius = 3;
    decomposition = 0;
    se = strel('disk', radius, decomposition);
    a = imopen(BW, se);
    
    % Threshold image - adaptive threshold
    %figure; imshow(a)
    
%     figure(80);
%     subplot(2, 2, [1 2]) % top row of the 3x3 grid
%     imagesc(a)
    stats = regionprops("table",a,"Centroid", ...
        "MajorAxisLength","MinorAxisLength");
    centers_t = stats.Centroid;
    diameters = mean([stats.MajorAxisLength stats.MinorAxisLength],2);
    radii_t = diameters/2;
%     hold on
%     viscircles(centers_t,radii_t)
%     plot(centers_t(:,1),centers_t(:,2),'+k') 
%     title('Detections using conventional method')

    % Detectamos en digital
    t_est_dig_copy=detrend_t_est_digcell{i};
    % figure; 
    % imshow(detrend_t_est_digcell{i},[])
    B = sort(detrend_t_est_digcell{i}(:));
    C = unique(B);
    detrend_t_est_digcell{i}(isnan(detrend_t_est_digcell{i}))=C(2); %mean(mean(t_sin_perfil,"omitnan"),"omitnan");
    detrend_t_est_digcell{i}(detrend_t_est_digcell{i}==0)=C(2); %mean(mean(t_sin_perfil,"omitnan"),"omitnan");
    
    X=detrend_t_est_digcell{i};
    Xmin = min(X(:));
    Xmax = max(X(:));
    if isequal(Xmax,Xmin)
        X = 0*X;
    else
        X = (X - Xmin) ./ (Xmax - Xmin);
    end
    
    % Threshold image - manual threshold
    BW = X > 4.588200e-01;
    
    % Open mask with disk
    radius = 3;
    decomposition = 0;
    se = strel('disk', radius, decomposition);
    a = imopen(BW, se);
    
    % Threshold image - adaptive threshold
    
%     subplot(122)
%     imagesc(a)
    stats = regionprops("table",a,"Centroid", ...
        "MajorAxisLength","MinorAxisLength");
    centers_dig = stats.Centroid;
    diameters = mean([stats.MajorAxisLength stats.MinorAxisLength],2);
    radii_dig = diameters/2;
%     hold on
%     viscircles(centers_dig,radii_dig)
%     plot(centers_dig(:,1),centers_dig(:,2),'+k')
%     title('Detections using digital features')
%     sgtitle(['Pose ' num2str(i)])



    %%%% Ordenamos
    
    % Crear conjuntos de puntos (ejemplo)
    A = centers_dig;
    B = centers_t;
    
    % Ordenar B para minimizar la distancia a A usando el método greedy
    center_t_ordenado = ordenar_puntos_greedy(A, B);
    
    if size(center_t_ordenado,1)>63
        [~, idx]=max(center_t_ordenado(:,2));
        center_t_ordenado(idx,:)=[];
        centers_dig(idx,:)=[];
%         figure;
%         imagesc(a)
%         hold on
%         plot(center_t_ordenado(:,1),center_t_ordenado(:,2),'xk')
%         size(center_t_ordenado)
    end
    % Visualizar la correspondencia antes y después
    % figure;
    % subplot(1,2,1);
    % scatter(A(:,1), A(:,2), 'ro'); hold on;
    % scatter(B(:,1), B(:,2), 'bx');
    % title('Antes de ordenar');
    % 
    % subplot(1,2,2);
    % scatter(A(:,1), A(:,2), 'ro'); hold on;
    % scatter(center_t_ordenado(:,1), center_t_ordenado(:,2), 'bx');
    % title('Después de ordenar');
    
    
    subplot(2,3,1); 
    imagesc(t_est_trad_copy)
    hold on
    for j=1:size(center_t_ordenado,1)
        plot(A(j,1),A(j,2),'+b')
        plot(center_t_ordenado(j,1),center_t_ordenado(j,2),'xk')
        %rms_t=norm(centers_t-centers_dig)
        %pause(0.5);
    end




    rms_t(i)=(sum(vecnorm(center_t_ordenado-centers_dig,2,2)))/size(center_t_ordenado,1)
    subplot(2,3,[2 3]); 
    plot(rms_t(1:i),'-ob')
    xlabel('Poses')
    ylabel('RMSE (px)')
    xlim([0 13])
    
    errormapt=NaN.*ones(size(t_est_trad_copy));
    error_pixel=vecnorm(center_t_ordenado-centers_dig,2,2);
    for k=1:size(center_t_ordenado,1)
        errormapt(floor(center_t_ordenado(k,2)),floor(center_t_ordenado(k,1)))=error_pixel(k);
    end
    %     figure; imagesc(errormapt)
    %     colorbar
    
    %     figure;
    %     stem(error_pixel,'ob')
    
    % dif_centers=centers_dig-center_t_ordenado;
    % figure;
    % scatter(dif_centers(:,1),dif_centers(:,2))
    
    
    % figure;
    % stem(center_t_ordenado(:,2))
    % hold on
    % stem(centers_dig(:,2))
    
    
    [x_error,y_error]=meshgrid(1:1280,1:1024);
    x_error=x_error(:);
    y_error=y_error(:);
    errormapt_interpol= griddata(double(center_t_ordenado(:,1)),double(center_t_ordenado(:,2)),error_pixel,double(x_error),double(y_error),"cubic");
    subplot(2,3,[4]);  
    imagesc(reshape(errormapt_interpol,[1024, 1280]))
    colorbar
    sgtitle(['Pose ' num2str(i)])
    pause(1.5);

end
pause();
subplot(2,3,[2 3]); 
plot(rms_t,'-ob')
xlabel('Poses')
ylabel('RMSE (px)')
xlim([0 13])
hold on
yl = yline(mean(rms_t),'--',['Mean: ' num2str(round(mean(rms_t),2)) ' px.'],'LabelVerticalAlignment','bottom','LineWidth',1);
yl.FontSize = 7;

%%
close all
figure;
subplot(121)
imagesc(detrend_z_est_digcell{1})
hold on
yline(300,'--')
xline(550,'--')

subplot(122)
imagesc(detrend_t_est_digcell{1})
hold on
yline(300,'--')
xline(550,'--')

c = colorbar; 
c.Label.String = '° C';


figure;
subplot(211)
yyaxis left
plot(detrend_z_est_digcell{1}(300,:),'b')
ylabel('Z (mm)')
hold on
plot(detrend_z_est_tradcell{1}(300,:),'-k')

yyaxis right
plot(detrend_t_est_digcell{1}(300,:),'r')
xlabel('X (mm)')
ylabel('t (°C)')
title('Profile in x')
hold on
plot(detrend_t_est_tradcell{1}(300,:),'g')
xlim([100 1200])

subplot(212)
yyaxis left
plot((detrend_z_est_digcell{1}(10:600,550)),'b')
ylabel('Z (mm)')
hold on
plot((detrend_z_est_tradcell{1}(10:600,550)),'-k')
%ylim([-0.4 0.4])

yyaxis right
plot(detrend_t_est_digcell{1}(10:600,550),'r')
xlabel('X (mm)')
ylabel('t (°C)')
title('Profile in y')
hold on
plot(detrend_t_est_tradcell{1}(10:600,550),'g')
xlim([0 600])


%% Perfil vertical
close all
figure;
imagesc(t_est_digcell{1})
hold on
plot(550,10:600,'.k')


figure;
plot(z_est_digcell{1}(10:600,550))

t_est=t_est_digcell{1}(10:600,550);
z_est=z_est_digcell{1}(10:600,550);
tam=length(t_est);
%tam=length(X_trad(n_profile,:));

% Vista de perfiles en la misma grafica sin perfil primario
profile_real=-z_est;
mask=isnan(profile_real);
profile=profile_real(~mask);



p = polyfit(1:length(profile),profile',2); 
f = polyval(p,1:length(profile)); 
profile_z_sin_primary=profile-f';
% figure;
% plot(f)
% hold on
% plot(profile)

%quitamos perfil primario
%Buscamos el perfil termico
profile_t=t_est;
profile_t=profile_t(~mask);

size(profile_z_sin_primary);
size(profile_t);
figure(7);
subplot(212)
yyaxis left
offset=-min(profile_z_sin_primary);
offset=0;
plot((profile_z_sin_primary+1*offset),'b')
ylabel('Z (mm)')
% yt = get(gca, 'YTick');
% set(gca, 'YTick',yt, 'ZTickLabel',fliplr(yt))
%ylim([-1 1])
yyaxis right
plot(profile_t)
xlabel('X (mm)')
ylabel('t (°C)')
title('Profile in x')



%quitamos perfil secundario
p1 = polyfit(1:length(profile_z_sin_primary),profile_z_sin_primary',3); 
f2 = polyval(p1,1:length(profile_z_sin_primary)); 
profile_z_final=profile_z_sin_primary-f2';
figure(8);
subplot(212)
yyaxis left
plot((profile_z_sin_primary-f2))
ylabel('Z (mm)')
%ylim([-1 1])
yyaxis right
plot(profile_t)
xlabel('X (mm)')
ylabel('t (°C)')
title('Traditional')
grid on
% if ~isempty(profile_z_sin_primary-f2)
%     z_est_final=NaN*ones(1,tam);
%     t_est_final=NaN*ones(1,tam);
%     z_est_final(~mask)=(profile_z_sin_primary-f2);
%     t_est_final(~mask)=profile_t;
% else
%     z_est_final=NaN*ones(1,tam);
%     t_est_final=NaN*ones(1,tam);
% end

if ~isempty(profile_z_sin_primary)
    z_est_final=NaN*ones(1,tam);
    t_est_final=NaN*ones(1,tam);
    z_est_final(~mask)=(profile_z_sin_primary);
    t_est_final(~mask)=profile_t;
else
    z_est_final=NaN*ones(1,tam);
    t_est_final=NaN*ones(1,tam);
end
figure;
yyaxis left
plot(z_est_final)
yyaxis right
plot(t_est_final)

figure(9);

aux=(profile_z_sin_primary);
tam_aux=1:length(aux);
%tam_aux=550:700;
%tam_aux=570:590
d31=aux(tam_aux);
t31=profile_t(tam_aux);

yyaxis left
plot(aux(tam_aux))
ylabel('Z (mm)')
%ylim([-1 1])
yyaxis right
plot(profile_t(tam_aux))
xlabel('X (px)')
ylabel('t (°C)')
%title('Digital')
grid on
hold on
%%%%%%%


figure;
plot(z_est_tradcell{1}(10:600,550))

t_est=t_est_tradcell{1}(10:600,550);
z_est=z_est_tradcell{1}(10:600,550);
tam=length(t_est);
%tam=length(X_trad(n_profile,:));

% Vista de perfiles en la misma grafica sin perfil primario
profile_real=-z_est;
mask=isnan(profile_real);
profile=profile_real(~mask);



p = polyfit(1:length(profile),profile',2); 
f = polyval(p,1:length(profile)); 
profile_z_sin_primary=profile-f';
% figure;
% plot(f)
% hold on
% plot(profile)

%quitamos perfil primario
%Buscamos el perfil termico
profile_t=t_est;
profile_t=profile_t(~mask);

size(profile_z_sin_primary);
size(profile_t);
figure(7);
subplot(212)
yyaxis left
offset=-min(profile_z_sin_primary);
offset=0;
plot((profile_z_sin_primary+1*offset),'b')
ylabel('Z (mm)')
% yt = get(gca, 'YTick');
% set(gca, 'YTick',yt, 'ZTickLabel',fliplr(yt))
%ylim([-1 1])
yyaxis right
plot(profile_t)
xlabel('X (mm)')
ylabel('t (°C)')
title('Profile in x')



%quitamos perfil secundario
p1 = polyfit(1:length(profile_z_sin_primary),profile_z_sin_primary',3); 
f2 = polyval(p1,1:length(profile_z_sin_primary)); 
profile_z_final=profile_z_sin_primary-f2';
figure(8);
subplot(212)
yyaxis left
plot((profile_z_sin_primary))
ylabel('Z (mm)')
%ylim([-1 1])
yyaxis right
plot(profile_t)
xlabel('X (mm)')
ylabel('t (°C)')
title('Traditional')
grid on
% if ~isempty(profile_z_sin_primary-f2)
%     z_est_final=NaN*ones(1,tam);
%     t_est_final=NaN*ones(1,tam);
%     z_est_final(~mask)=(profile_z_sin_primary-f2);
%     t_est_final(~mask)=profile_t;
% else
%     z_est_final=NaN*ones(1,tam);
%     t_est_final=NaN*ones(1,tam);
% end

if ~isempty(profile_z_sin_primary)
    z_est_final=NaN*ones(1,tam);
    t_est_final=NaN*ones(1,tam);
    z_est_final(~mask)=(profile_z_sin_primary);
    t_est_final(~mask)=profile_t;
else
    z_est_final=NaN*ones(1,tam);
    t_est_final=NaN*ones(1,tam);
end
figure;
yyaxis left
plot(z_est_final)
yyaxis right
plot(t_est_final)

figure(9);
aux=(profile_z_sin_primary);


d32=aux(tam_aux);
t32=profile_t(tam_aux);

yyaxis left
plot(aux(tam_aux),'-k')
ylabel('Z (mm)')
%ylim([-1 1])
yyaxis right
plot(profile_t(tam_aux),'-g')
xlabel('X (px)')
ylabel('t (°C)')
%title('Digital')
grid on
legend('Z  by digital features','Z  by conventional approach','t by digital features','t by conventional approach')
grid minor

%%
size(t31)
size(t32)
rmse=sqrt(sum((t31-t32).^2)/size(t32,1))

figure; plot(abs(t31-t32))
hold on
plot(t31)

%%
figure;
subplot(141)
name_texture=[direc 'pose_0' num2str(poses(1)) '\im_51.png'];
texture=imread(name_texture);
imshow(texture)

subplot(141)
imagesc(i)

%% Perfil Horizontal
close all
figure;
imagesc(t_est_digcell{1})
hold on
plot(100:1280,310,'.k')
figure;
plot(z_est_digcell{1}(310,100:end))

t_est=t_est_digcell{1}(310,100:end);
z_est=z_est_digcell{1}(310,100:end);
tam=length(t_est);
%tam=length(X_trad(n_profile,:));

% Vista de perfiles en la misma grafica sin perfil primario
profile_real=-z_est;
mask=isnan(profile_real);
profile=profile_real(~mask);



p = polyfit(1:length(profile),profile,2); 
f = polyval(p,1:length(profile)); 
profile_z_sin_primary=profile-f;
% figure;
% plot(f)
% hold on
% plot(profile)

%quitamos perfil primario
%Buscamos el perfil termico
profile_t=t_est;
profile_t=profile_t(~mask);

size(profile_z_sin_primary);
size(profile_t);
figure(7);
subplot(212)
yyaxis left
offset=-min(profile_z_sin_primary);
offset=0;
plot((profile_z_sin_primary+1*offset),'b')
ylabel('Z (mm)')
% yt = get(gca, 'YTick');
% set(gca, 'YTick',yt, 'ZTickLabel',fliplr(yt))
%ylim([-1 1])
yyaxis right
plot(profile_t)
xlabel('X (mm)')
ylabel('t (°C)')
title('Profile in x')



%quitamos perfil secundario
p1 = polyfit(1:length(profile_z_sin_primary),profile_z_sin_primary,2); 
f2 = polyval(p1,1:length(profile_z_sin_primary)); 
profile_z_final=profile_z_sin_primary-f2;
figure(8);
subplot(212)
yyaxis left
plot((profile_z_sin_primary-f2))
ylabel('Z (mm)')
%ylim([-1 1])
yyaxis right
plot(profile_t)
xlabel('X (mm)')
ylabel('t (°C)')
title('Traditional')
grid on
% if ~isempty(profile_z_sin_primary-f2)
%     z_est_final=NaN*ones(1,tam);
%     t_est_final=NaN*ones(1,tam);
%     z_est_final(~mask)=(profile_z_sin_primary-f2);
%     t_est_final(~mask)=profile_t;
% else
%     z_est_final=NaN*ones(1,tam);
%     t_est_final=NaN*ones(1,tam);
% end

if ~isempty(profile_z_sin_primary)
    z_est_final=NaN*ones(1,tam);
    t_est_final=NaN*ones(1,tam);
    z_est_final(~mask)=(profile_z_sin_primary);
    t_est_final(~mask)=profile_t;
else
    z_est_final=NaN*ones(1,tam);
    t_est_final=NaN*ones(1,tam);
end
figure;
yyaxis left
plot(z_est_final)
yyaxis right
plot(t_est_final)

figure(10);

aux=(profile_z_sin_primary);
tam_aux=1:length(aux);
%tam_aux=550:700;
%tam_aux=570:590


yyaxis left
plot(aux(tam_aux))
ylabel('Z (mm)')
%ylim([-1 1])
yyaxis right
plot(profile_t(tam_aux))
xlabel('X (px)')
ylabel('t (°C)')
%title('Digital')
grid on
hold on
%%%%%%%


figure;
plot(z_est_tradcell{1}(310,100:end))

t_est=t_est_tradcell{1}(310,100:end);
z_est=z_est_tradcell{1}(310,100:end);
tam=length(t_est);
%tam=length(X_trad(n_profile,:));

% Vista de perfiles en la misma grafica sin perfil primario
profile_real=-z_est;
mask=isnan(profile_real);
profile=profile_real(~mask);



p = polyfit(1:length(profile),profile,2); 
f = polyval(p,1:length(profile)); 
profile_z_sin_primary=profile-f;
% figure;
% plot(f)
% hold on
% plot(profile)

%quitamos perfil primario
%Buscamos el perfil termico
profile_t=t_est;
profile_t=profile_t(~mask);

size(profile_z_sin_primary);
size(profile_t);
figure(7);
subplot(212)
yyaxis left
offset=-min(profile_z_sin_primary);
offset=0;
plot((profile_z_sin_primary+1*offset),'b')
ylabel('Z (mm)')
% yt = get(gca, 'YTick');
% set(gca, 'YTick',yt, 'ZTickLabel',fliplr(yt))
%ylim([-1 1])
yyaxis right
plot(profile_t)
xlabel('X (mm)')
ylabel('t (°C)')
title('Profile in x')



%quitamos perfil secundario
p1 = polyfit(1:length(profile_z_sin_primary),profile_z_sin_primary,2); 
f2 = polyval(p1,1:length(profile_z_sin_primary)); 
profile_z_final=profile_z_sin_primary-f2;
figure(8);
subplot(212)
yyaxis left
plot((profile_z_sin_primary))
ylabel('Z (mm)')
%ylim([-1 1])
yyaxis right
plot(profile_t)
xlabel('X (mm)')
ylabel('t (°C)')
title('Traditional')
grid on
% if ~isempty(profile_z_sin_primary-f2)
%     z_est_final=NaN*ones(1,tam);
%     t_est_final=NaN*ones(1,tam);
%     z_est_final(~mask)=(profile_z_sin_primary-f2);
%     t_est_final(~mask)=profile_t;
% else
%     z_est_final=NaN*ones(1,tam);
%     t_est_final=NaN*ones(1,tam);
% end

if ~isempty(profile_z_sin_primary)
    z_est_final=NaN*ones(1,tam);
    t_est_final=NaN*ones(1,tam);
    z_est_final(~mask)=(profile_z_sin_primary);
    t_est_final(~mask)=profile_t;
else
    z_est_final=NaN*ones(1,tam);
    t_est_final=NaN*ones(1,tam);
end
figure;
yyaxis left
plot(z_est_final)
yyaxis right
plot(t_est_final)

figure(10);
aux=(profile_z_sin_primary);

yyaxis left
plot(aux(tam_aux),'-k')
ylabel('Z (mm)')
%ylim([-1 1])
yyaxis right
plot(profile_t(tam_aux),'-g')
xlabel('X (px)')
ylabel('t (°C)')
%title('Digital')
grid on
legend('Z  by Digital features','Z  by Traditional','t by Digital features','t by Traditional')
grid minor
%% Imagen cualitativa
close all
clc
%clear all
direc='D:\R3D\Adquisiciones\Experimento 54\Objects\SL\';
aux=num2str(0);
% Cargamos digital features image
pc_dig = pcread([direc 'Digital Features 1 UTB\new_obj_vis_ir_0' aux '.ply']);
figure;
%subplot(131)
pcshow(pc_dig,'BackgroundColor','w');
xlabel('x (mm)')
ylabel('y (mm)')
zlabel('z (mm)')
c = colorbar; 
c.Label.String = '° C';
set(gca, 'Zdir', 'reverse')
set(gca, 'Xdir', 'reverse')
view(30,45)
grid off
zlim([350 430])
%title('Digital features I')

% % Cargamos traditional image
% aux=num2str(1);
% 
% pc_trad = pcread([direc 'Digital Features 1 UTB\new_obj_vis_ir_0' aux '.ply']);
% subplot(132)
% pcshow(pc_trad,'BackgroundColor','w');
% xlabel('x (mm)')
% ylabel('y (mm)')
% zlabel('z (mm)')
% c = colorbar; 
% c.Label.String = '° C';
% set(gca, 'Zdir', 'reverse')
% set(gca, 'Xdir', 'reverse')
% view(30,45)
% 
% grid off
% zlim([350 430])
% title('Traditional')
% 
% 
% 
% % Cargamos traditional image
% aux=num2str(0);
% direc='D:\R3D\Adquisiciones\Experimento 53\Objects\SL\';
% pc_trad = pcread([direc 'Digital Features 1 UTB\new_obj_vis_ir_0' aux '.ply']);
% subplot(133)
% pcshow(pc_trad,'BackgroundColor','w');
% xlabel('x (mm)')
% ylabel('y (mm)')
% zlabel('z (mm)')
% c = colorbar; 
% c.Label.String = '° C';
% set(gca, 'Zdir', 'reverse')
% set(gca, 'Xdir', 'reverse')
% view(30,45)
% 
% grid off
% zlim([350 430])
% title('Traditional')

figure; %subplot(131)
im0=imread('C:\Users\Eberto Benjumea\Desktop\OneDrive - Universidad Tecnológica de Bolívar\Images for Events and Papers\Applied Optics\new_figures\54\ir Image.png');
imshow(im0)
% subplot(132)
% im1=imread('C:\Users\Eberto Benjumea\Desktop\OneDrive - Universidad Tecnológica de Bolívar\Images for Events and Papers\Applied Optics\new_figures\54\ir Image2.png');
% imshow(im1)
% subplot(133)
% im2=imread('C:\Users\Eberto Benjumea\Desktop\OneDrive - Universidad Tecnológica de Bolívar\Images for Events and Papers\Applied Optics\new_figures\53\00.png');
% imshow(im2)


% Digital features I
x_dig=pc_dig.Location(:,1);
y_dig=pc_dig.Location(:,2);
z_dig=pc_dig.Location(:,3);
t_dig=pc_dig.Intensity;
[X_dig,Y_dig]=meshgrid([min(x_dig):0.05:max(x_dig)], [min(y_dig):0.05:max(y_dig)]);

%     [X_dig,Y_dig]=meshgrid(linspace(min(x_dig),max(x_dig),1280), linspace(min(y_dig),max(y_dig),1024));



z_est_dig=griddata(x_dig,y_dig,z_dig,X_dig,Y_dig,'cubic');
t_est_dig=griddata(x_dig,y_dig,t_dig,X_dig,Y_dig,'cubic');
z_est_dig(z_est_dig>max(z_dig))=NaN;
t_est_dig(t_est_dig>max(t_dig))=NaN;
z_est_dig(z_est_dig<min(z_dig))=NaN;
t_est_dig(t_est_dig<min(t_dig))=NaN;



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



load([direc 'Digital Features 1 UTB\obj_vis_ir_0' aux '.mat'])

%% analisis de t respecto a Z linea
close all
for i=1:length(detrend_z_est_digcell)
    a=detrend_z_est_digcell{i};

    filename = sprintf('dig_%02d.mat', i); % Nombre del archivo
    mat = detrend_z_est_digcell{i}; % Extraer matriz
    save([direc filename], 'mat'); % Guardar en archivo .mat

end

close all
for i=1:length(detrend_z_est_tradcell)
    a=detrend_z_est_tradcell{i};

    filename = sprintf('trad_%02d.mat', i); % Nombre del archivo
    mat = detrend_z_est_tradcell{i}; % Extraer matriz
    save([direc filename], 'mat'); % Guardar en archivo .mat

end


for i=1:length(detrend_z_est_digcell)
    a=detrend_z_est_digcell{i};

    filename = sprintf('dig_%02d.png', i); % Nombre del archivo
    img = mat2gray(a); % Normalizar valores entre 0 y 1 (importante para imágenes)
    imwrite(img, [direc filename]); % Guardar como imagen PNG




end

close all
for i=1:length(detrend_z_est_tradcell)
    a=detrend_z_est_tradcell{i};

    filename = sprintf('trad_%02d.png', i); % Nombre del archivo
    img = mat2gray(a); % Normalizar valores entre 0 y 1 (importante para imágenes)
    imwrite(img, [direc filename]); % Guardar como imagen PNG

end

    filename = sprintf('matrix_%02d.png', i); % Nombre del archivo
    img = mat2gray(C{i}); % Normalizar valores entre 0 y 1 (importante para imágenes)
    imwrite(img, filename); % Guardar como imagen PNG




a=detrend_z_est_digcell{1};
mean_dig=mean(mean(a,'omitnan'),'omitnan');
a(isnan(a))=mean_dig;
b=edge(a,"log",0.0025);
figure; imshow(b)


[b, threshOut]=edge(a,"canny");
figure; imshow(b)
threshOut
close all
[b, threshOut]=edge(a,"canny",[0.0125    0.025]);
figure; imshow(b)

[B, L] = bwboundaries(b, 'noholes'); % Encuentra los contornos en edge_image
figure;
imshow(L)
hold on
for k = 1:length(B)
   boundary = B{k};
   plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2)
end
% Threshold image - adaptive threshold
figure; imshow(a)
[mask_w_w,numWhite] = bwlabel(a);
mask_w_w = bwareaopen(mask_w_w,3);  %Retirando manchas blancas
figure; imshow(mask_w_w)

%% Analisis de t respecto a Z

% Plano ideal
close all
figure;
detrend_z_est_digcell={};
detrend_t_est_digcell={};
for i=1:length(z_est_digcell)
    imagesc(z_est_digcell{i})
        mean_dig=mean(mean(z_est_digcell{i},'omitnan'),'omitnan');
        aux_dig=z_est_digcell{i}-mean_dig;
        aux_digx=X_digcell{i};
        aux_digy=Y_digcell{i};

        x_aux=aux_digx(:);
        y_aux=aux_digy(:);
        z_aux=aux_dig(:);
        x=x_aux(~isnan(z_aux));
        y=y_aux(~isnan(z_aux));
        z=z_aux(~isnan(z_aux));
        DM = [x, y, ones(size(z))];                             % Design Matrix
        %B = DM\z;                                               % Estimate Parameters
        B=lsqminnorm(DM,z);
        Z = B(1)*x + B(2)*y + B(3)*ones(size(x));


%         figure;
%         plot3(x,y,Z,'.r')
%         view([45 45])
%         title('Reconstructed ideal 3D plane')
%         xlabel('x (mm)')
%         ylabel('y (mm)')
%         zlabel('z (mm)')
% 
%         figure;
%         plot3(x,y,z-Z,'.r')
%         view([45 45])
%         title('Reconstructed ideal 3D plane')
%         xlabel('x (mm)')
%         ylabel('y (mm)')
%         zlabel('z (mm)')

        [X,Y] = meshgrid(linspace(min(x),max(x),1280), linspace(min(y),max(y),1024));


            [X_trad,Y_trad]=meshgrid(linspace(min(x_trad),max(x_trad),1280), linspace(min(y_trad),max(y_trad),1024));
    
    
    
        z_sin_media=griddata(x,y,z-Z,X,Y,'cubic');
        s = surf(-X_trad,Y_trad,-z_sin_media,'FaceColor', 'r',...
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
        mean_dig=mean(mean(z_sin_media,'omitnan'),'omitnan');
        z_sin_media(isnan(z_sin_media))=mean_dig;
        detrend_z_est_digcell{i}=z_sin_media;


%         [z_est_final,t_est_final]=detrend_profile(aux_dig, t_est_digcell{i}(j,:));
%         detrend_z_est_digcell{i}(j,:)=z_est_final;
%         detrend_t_est_digcell{i}(j,:)=t_est_final;        


    imagesc(detrend_z_est_digcell{i})
    title(['Pose ' num2str(i)])
    pause(1);
end


a=detrend_z_est_digcell{1};
a=imadjust(a);
figure; imagesc(a)
colorbar

b=edge(a,"log");
figure; imshow(b)


b=edge(a,"canny");
figure; imshow(b)



% Threshold image - adaptive threshold
figure; imshow(a)
[mask_w_w,numWhite] = bwlabel(a);
mask_w_w = bwareaopen(mask_w_w,3);  %Retirando manchas blancas
figure; imshow(mask_w_w)





detrend_z_est_tradcell={};
detrend_t_est_tradcell={};
for i=1:length(z_est_tradcell)
    imagesc(z_est_tradcell{i})
        mean_trad=mean(mean(z_est_tradcell{i},'omitnan'),'omitnan');
        aux_trad=z_est_tradcell{i}-mean_trad;
        aux_tradx=X_tradcell{i};
        aux_trady=Y_tradcell{i};

        x_aux=aux_tradx(:);
        y_aux=aux_trady(:);
        z_aux=aux_trad(:);
        x=x_aux(~isnan(z_aux));
        y=y_aux(~isnan(z_aux));
        z=z_aux(~isnan(z_aux));
        DM = [x, y, ones(size(z))];                             % Design Matrix
        %B = DM\z;                                               % Estimate Parameters
        B=lsqminnorm(DM,z);
        Z = B(1)*x + B(2)*y + B(3)*ones(size(x));


%         figure;
%         plot3(x,y,Z,'.r')
%         view([45 45])
%         title('Reconstructed ideal 3D plane')
%         xlabel('x (mm)')
%         ylabel('y (mm)')
%         zlabel('z (mm)')
% 
%         figure;
%         plot3(x,y,z-Z,'.r')
%         view([45 45])
%         title('Reconstructed ideal 3D plane')
%         xlabel('x (mm)')
%         ylabel('y (mm)')
%         zlabel('z (mm)')

        [X,Y] = meshgrid(linspace(min(x),max(x),1280), linspace(min(y),max(y),1024));


            [X_trad,Y_trad]=meshgrid(linspace(min(x_trad),max(x_trad),1280), linspace(min(y_trad),max(y_trad),1024));
    
    
    
        z_sin_media=griddata(x,y,z-Z,X,Y,'cubic');
%         s = surf(-X_trad,Y_trad,-z_sin_media,'FaceColor', 'r',...
%                     'EdgeColor', 'none',...
%                     'FaceLighting', 'phong');
%         
%         
%         set(gca, 'DataAspectRatio', [1, 1, 1])
%         axis equal;
%         view(180,90);
%         camlight right
%         axis on
%         grid on
%         xlabel('x (mm)')
%         ylabel('y (mm)')
%         zlabel('z (mm)')
        detrend_z_est_tradcell{i}=z_sin_media;


%         [z_est_final,t_est_final]=detrend_profile(aux_trad, t_est_tradcell{i}(j,:));
%         detrend_z_est_tradcell{i}(j,:)=z_est_final;
%         detrend_t_est_tradcell{i}(j,:)=t_est_final;        


    imagesc(detrend_z_est_tradcell{i})
    title(['Pose ' num2str(i)])
    pause(0.5);
end

figure;
yyaxis left
plot(detrend_z_est_tradcell{1}(300,:),'b')
ylabel('Z (mm)')
hold on
plot(detrend_z_est_digcell{1}(300,:),'-k')




%%
close all
poses=[5 8 13 14 16 17 18 19 20 21 24 25 26 28 29 30 31 32 33 35 36 39 40 41 42 43 44 46 47 49 53 60];
length(poses)
i=1
    if poses(i)<10
        aux=['0' num2str(poses(i))];
    else
        aux=num2str(poses(i));
    end
    disp(['Pose ' aux])
    figure;

    % Cargamos digital features image
    pc_dig = pcread([direc 'Digital Features 1 UTB\new_obj_vis_ir_' aux '.ply']);
    pcshow(pc_dig);
    xlabel('x (mm)')
    ylabel('y (mm)')
    zlabel('z (mm)')
    c = colorbar; 
    c.Label.String = '° C';
    set(gca, 'Zdir', 'reverse')
    set(gca, 'Xdir', 'reverse')
    view(30,45)
    grid off
    zlim([350 430])
    title('Digital features I')
    
    % Cargamos traditional image
 
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
%     


