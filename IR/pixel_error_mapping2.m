close all
clc
clear all
direc='D:\R3D\Adquisiciones\Experimento 54\Map_cal\SL\';
n_poses=62;
addpath 'C:\Users\Eberto Benjumea\Desktop\OneDrive - Universidad Tecnológica de Bolívar\Doctorado\Proyectos\Codigos Reconstruccion 3D Learning\Codigos Tesis All Digital Features'


%%
close all
%poses=[5 8 18 19 20 24 28 29 30 33 41 42 47];
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



end

figure;
for i=1:length(z_est_digcell)
    if poses(i)<10
        aux=['0' num2str(poses(i))];
    else
        aux=num2str(poses(i));
    end

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




end












%%
for i=2:size(detrend_z_est_digcell,2)  % 2,3,7,8
    i=13
    close all
     t_est_trad=detrend_t_est_tradcell{i};
    t_est_trad_copy=t_est_trad;
    figure; 
    imshow(t_est_trad,[])
    B = sort(t_est_trad(:));
    C = unique(B);
    t_est_trad(isnan(t_est_trad))=C(2); %mean(mean(t_sin_perfil,"omitnan"),"omitnan");
    t_est_trad(t_est_trad==0)=C(2); %mean(mean(t_sin_perfil,"omitnan"),"omitnan");
    
    X=t_est_trad;
    Xmin = min(X(:));
    Xmax = max(X(:));
    if isequal(Xmax,Xmin)
        X = 0*X;
    else
        X = (X - Xmin) ./ (Xmax - Xmin);
    end
    
    % Threshold image - adaptive threshold
    BW = imbinarize(X, 'adaptive', 'Sensitivity', 0.100000, 'ForegroundPolarity', 'bright');
    figure; imshow(BW)
    % Open mask with disk
    radius = 10;
    decomposition = 0;
    se = strel('disk', radius, decomposition);
    a = imopen(BW, se);
    
    % Threshold image - adaptive threshold
    figure; imshow(a)
    [mask_w_w,numWhite] = bwlabel(a);
    mask_w_w = bwareaopen(mask_w_w,2500);  %Retirando manchas blancas
    figure; imshow(mask_w_w)
    
    %a=(imbinarize(t_sin_perfil));
    a=mask_w_w ;
    b= ~bwareaopen(a,25000);   %Retirando manchas negras
    a=logical(a.*b);
     figure; 
     imshow(a,[])

    
    figure;
    imagesc(a)
    stats = regionprops("table",a,"Centroid", ...
        "MajorAxisLength","MinorAxisLength","Circularity")
    
    center_aux=stats.Centroid(stats.Circularity>0.8,:);
%     hold on
%     plot(center_aux(:,1),center_aux(:,2),'+r')
    stats(stats.Circularity<0.8,:)=[];
    diameters = mean([stats.MajorAxisLength stats.MinorAxisLength],2);
    radii_t = diameters/2;
    centers_t_trad = stats.Centroid;
    hold on
    viscircles(centers_t_trad,radii_t)
    plot(centers_t_trad(:,1),centers_t_trad(:,2),'+k')

    %

      t_est_dig=detrend_t_est_digcell{i};
    t_est_dig_copy=t_est_dig;
    figure; 
    imshow(t_est_dig,[])
    B = sort(t_est_dig(:));
    C = unique(B);
    t_est_dig(isnan(t_est_dig))=C(2); %mean(mean(t_sin_perfil,"omitnan"),"omitnan");
    t_est_dig(t_est_dig==0)=C(2); %mean(mean(t_sin_perfil,"omitnan"),"omitnan");
    
    X=t_est_dig;
    Xmin = min(X(:));
    Xmax = max(X(:));
    if isequal(Xmax,Xmin)
        X = 0*X;
    else
        X = (X - Xmin) ./ (Xmax - Xmin);
    end
    
    % Threshold image - adaptive threshold
    BW = imbinarize(X, 'adaptive', 'Sensitivity', 0.100000, 'ForegroundPolarity', 'bright');
    figure; imshow(BW)
    % Open mask with disk
    radius = 10;
    decomposition = 0;
    se = strel('disk', radius, decomposition);
    a = imopen(BW, se);
    
    % Threshold image - adaptive threshold
    figure; imshow(a)
    [mask_w_w,numWhite] = bwlabel(a);
    mask_w_w = bwareaopen(mask_w_w,2500);  %Retirando manchas blancas
    figure; imshow(mask_w_w)
    
    %a=(imbinarize(t_sin_perfil));
    a=mask_w_w ;
    b= ~bwareaopen(a,25000);   %Retirando manchas negras
    a=logical(a.*b);
     figure; 
     imshow(a,[])

    
    figure;
    imagesc(a)
    stats = regionprops("table",a,"Centroid", ...
        "MajorAxisLength","MinorAxisLength","Circularity")
    
    center_aux=stats.Centroid(stats.Circularity>0.8,:);
%     hold on
%     plot(center_aux(:,1),center_aux(:,2),'+r')
    stats(stats.Circularity<0.8,:)=[];
    diameters = mean([stats.MajorAxisLength stats.MinorAxisLength],2);
    radii_t = diameters/2;
    centers_t_dig = stats.Centroid;
    hold on
    viscircles(centers_t_dig,radii_t)
    plot(centers_t_dig(:,1),centers_t_dig(:,2),'+k')






    figure;
    imagesc(detrend_z_est_digcell{i})
    aux_z_orig=detrend_z_est_digcell{i};
    bCanny = edge(aux_z_orig, 'canny');  % Umbrales automáticos
    figure;
    imshow(bCanny)
    [centers_z_dig, radii, metric] = imfindcircles(bCanny,[25 50],'ObjectPolarity','dark',Sensitivity=0.928);
    figure;
    imagesc(bCanny)
    hold on
    plot(centers_z_dig(:,1),centers_z_dig(:,2),'ok')
    viscircles(centers_z_dig,radii)

    
    P1=centers_z_dig;
    P2=centers_t_dig;
    [idx, dist]  = knnsearch(P2, P1);          % vecino más próximo
    P2_sorted    = P2(idx,:);
    
    %--- visualizar ---
    figure;
    scatter(P1(:,1), P1(:,2), 60, 'b', 'filled'); hold on
    scatter(P2(:,1), P2(:,2), 60, 'r');         % sin ordenar
    quiver(P1(:,1), P1(:,2), ...
           P2_sorted(:,1)-P1(:,1), P2_sorted(:,2)-P1(:,2), ...
           0, 'k');                             % flechas de correspondencia
    legend({'P1','P2 (raw)','match'}, 'Location','best');
    title('Emparejamiento por vecino más cercano');
    axis equal
    diff_pixel=vecnorm(centers_z_dig-P2_sorted,2,2);
    figure; plot(diff_pixel)
    figure; histogram(diff_pixel)
    mean(diff_pixel)
    %
    figure;
    imagesc(detrend_z_est_tradcell{i})
    aux_z_orig=detrend_z_est_tradcell{i};
    bCanny = edge(aux_z_orig, 'canny');  % Umbrales automáticos
    figure;
    imagesc(bCanny)
    [centers_z_trad, radii, metric] = imfindcircles(bCanny,[25 50],'ObjectPolarity','dark',Sensitivity=0.928);
    figure;
    imagesc(aux_z_orig)
    hold on
    plot(centers_z_trad(:,1),centers_z_trad(:,2),'ok')
    viscircles(centers_z_trad,radii)
    P1=centers_z_trad;
    P2=centers_t_trad;
    [idx, dist]  = knnsearch(P2, P1);          % vecino más próximo
    P2_sorted    = P2(idx,:);
    
    %--- visualizar ---
    figure;
    scatter(P1(:,1), P1(:,2), 60, 'b', 'filled'); hold on
    scatter(P2(:,1), P2(:,2), 60, 'r');         % sin ordenar
    quiver(P1(:,1), P1(:,2), ...
           P2_sorted(:,1)-P1(:,1), P2_sorted(:,2)-P1(:,2), ...
           0, 'k');                             % flechas de correspondencia
    legend({'P1','P2 (raw)','match'}, 'Location','best');
    title('Emparejamiento por vecino más cercano');
    axis equal
    diff_pixel=vecnorm(centers_z_trad-P2_sorted,2,2);
    figure; plot(diff_pixel)
    figure; histogram(diff_pixel)
    mean(diff_pixel)





    pause();

end
