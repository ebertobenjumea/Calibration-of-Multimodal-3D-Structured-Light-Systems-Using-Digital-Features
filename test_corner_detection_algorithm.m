close all
clc
clear_ask=input('Do you want to clear the previously processed data? Yes:1  No:2\nYour answer:');
if clear_ask==1
    clear all
    disp('All data was cleared')
else
    disp('You are working with previously processed data')
end
inst_ask=input('UTB or CITEDI data? UTB:1  CITEDI:2\nYour answer:');
if inst_ask==1
    disp('You are working with UTB data!')
    direc='D:\R3D\Adquisiciones\Experimento 54\Digital features\';
    name_texture='P18_F26';
elseif inst_ask==2
    disp('You are working with CITEDI data!')
    direc='D:\R3D\2023\multimodal\Experimento 44\Digital Features 2\';
    %direc='D:\R3D\2023\multimodal\Experimento 18\';
    name_texture='P18_F26';
end

% Input data
n_poses=81;
n_corners=6;
camera1='camera 1\';
camera2='camera 2\';
IR='IR\';
format='.PNG';
cd('C:\Users\Eberto Benjumea\Desktop\OneDrive - Universidad Tecnológica de Bolívar\Doctorado\Proyectos\Codigos Reconstruccion 3D Learning\Codigos Tesis All Digital Features')
%% Mirror corner detection and mask creation for C1
clc
disp('Mirror corner detection and mask creation.')
disp(' ')


figure;
texture_img={};
mask={};
mask_final={};
corner_store={};
disp('Pose: ')
for i=1:n_poses
    if i<11
        aux=['0' num2str(i-1)];
    else
        aux=num2str(i-1);
    end
    disp(aux)
    % Leemos la imagen
    texture_img{i}=imread([direc camera1 'pose_' aux '\H\' name_texture format]);    %Accedemos a imagen de textura
    img=texture_img{i};
    width=size(img,2);
    height=size(img,1);
    % Buscamos las esquinas y la mascara via umbralziación
    sgtitle(['Pose ' num2str(i-1)])
    [mask{i},corner_store{i}]=find_mirror_nmask(img, n_corners,1);
    pause(0.5)
    
end
close all


%%
figure;
for i=1:n_poses
    
    if i<11
        aux=['0' num2str(i-1)];
    else
        aux=num2str(i-1);
    end
    disp(aux)
    % Leemos la imagen
    texture_img{i}=imread([direc camera1 'pose_' aux '\H\' name_texture format]);    %Accedemos a imagen de textura
    img=texture_img{i};




    % corner_contour_store=[];
    % corner_store=[];
    % c_semifinal_corners=[];
    % final_corners=[];
    img= imadjust(img);
    % Buscamos las esquinas y la mascara via umbralziación
    img=imsharpen(img);
    blur = imgaussfilt(img, 1); % '1' es el sigma, que puedes ajustar si es necesario
    level = graythresh(blur);        %Umbralizamos
    mask1 = imbinarize(blur,level);   
    [mask_w_w,numWhite] = bwlabel(mask1);
    mask_w_w = bwareaopen(mask_w_w,25000);  %Retirando manchas blancas
    [L,numBlack] = bwlabel(~mask_w_w); %Etiquetamos los agrupamientos de pixeles.
    mask= ~bwareaopen(L,25000);   %Retirando manchas negras
    
    % Paso 1: Encontrar contornos en la imagen binaria
    [B, L] = bwboundaries(mask, 'noholes'); % Encuentra los contornos en edge_image
    
    % Paso 2: Ordenar los contornos por área
    contour_areas = cellfun(@(x) polyarea(x(:,2), x(:,1)), B);
    [areas, sorted_indices] = sort(contour_areas, 'descend');
    sorted_contours = B(sorted_indices);
    
    % Paso 3: Iterar sobre los contornos ordenados para buscar el de 6 esquinas
    contour = sorted_contours{1};
    % Calcular el perímetro del contorno
    perimeter = sum(sqrt(sum(diff(contour).^2, 2)));
    % Aproximación del contorno utilizando el umbral de aproximación
    epsilon=2*perimeter/areas(1);
    if epsilon>0.2
        epsilon=0.02;
    end
    approx_contour = reducepoly(contour, epsilon);
    [C,ia,ic] = unique(approx_contour(:,1:2),'rows');
    approx_contour = approx_contour(ia,:);
    %disp(['Corners: ' num2str(size(approx_contour,1))])
    
    corner_contour_store=[approx_contour(:,2) approx_contour(:,1)];
    
    mirror_corners = pgonCorners(mask,n_corners,1000);
    corner_store=[mirror_corners(:,2) mirror_corners(:,1)];
    
    
    corners= detectHarrisFeatures(mask);
    corndata=corners.selectStrongest(500);
    corn_x=corndata.Location(:,1);
    corn_y=corndata.Location(:,2);
    posible_corners=[corn_x corn_y];
    
    cnt_valid_corner=0;
    s = regionprops(mask,'centroid'); 
    centroids = cat(1,s.Centroid);
    
    for k=1:size(corner_store,1)
        % Buscamos en corners Harris
        distancia=vecnorm(posible_corners-corner_store(k,:),2,2);
        [~, idx]=min(distancia);
        %harris_corners=[posible_corners(idx,:); corner_store{i}(k,:)];
        harris_corners=posible_corners(idx,:);
        % Buscamos en corners por contornos
        distancia=vecnorm(corner_contour_store-corner_store(k,:),2,2);
        [~, idx]=min(distancia);
        contour_corners=corner_contour_store(idx,:);
        semifinal_corners=[corner_store(k,:); contour_corners; harris_corners];
         
        dist_semifinal_corners=vecnorm(semifinal_corners-centroids,2,2);
        [maximo idx_few]=max(dist_semifinal_corners);
        final_corners(k,:)=semifinal_corners(idx_few,:);
    end
    disp(['Corners: ' num2str(size(final_corners,1))])


    imshow(img)
    hold on
    plot(contour(:,2),contour(:,1),'.b');







%     imshow(img)
%     hold on
%     plot(corner_store{i}(:,1),corner_store{i}(:,2),'+b');
%     subplot(122) 
%     imshow(mask,[])
    pause();
end
%% Busqueda de mascara en C2
figure;
texture_img2={};
mask2={};
mask_final2={};
corner_store2={};
disp('Pose: ')
for i=1:n_poses
    if i<11
        aux=['0' num2str(i-1)];
    else
        aux=num2str(i-1);
    end
    disp(aux)
    texture_img2{i}=imread([direc camera2 'pose_' aux '\H\' name_texture format]);    %Accedemos a imagen de textura
    img2=texture_img2{i};
    % Buscamos las esquinas y la mascara via umbralziación
    sgtitle(['Pose ' num2str(i-1)])
    [mask2{i},corner_store2{i}]=find_mirror_nmask(img2, n_corners,1);
    pause(0.5);
end
close all

%% Muestreo de C1 y C2
clc
close all

grid_size1=20;
[screen_points_c1,mask_points]=points_camera(texture_img{1},grid_size1);
grid_size2=3;
[screen_points_c2,mask_points2]=points_camera(texture_img2{1},grid_size2);

figure; 
subplot(121)
imshow(texture_img{1})
hold on
plot(screen_points_c1(:,1),screen_points_c1(:,2),'*r')
title(['Camera 1 samples with grid size= ' num2str(grid_size1)])
subplot(122) 
imshow(texture_img2{1})
hold on
plot(screen_points_c2(:,1),screen_points_c2(:,2),'.r')
title(['Camera 2 samples with grid size= ' num2str(grid_size2)])
pause(3);

%% Combinación de mascarás (Rectangulo y muestreo de C1)
close all
mask_ROI={};
mask_final=mask;
feature_points_c1_mirror={};
figure;
for i=1:length(mask_final)
    mask_ROI{i}=mask_final{i}.*mask_points;
    [ycoord, xcoord] = find(mask_ROI{i});
    feature_points_c1_mirror{i}=[xcoord, ycoord];
    subplot(121)
    imshow(texture_img{i})
    hold on
    plot(feature_points_c1_mirror{i}(:,1),feature_points_c1_mirror{i}(:,2),'.r')
    subplot(122)
    imshow(mask_final{i})
    hold on
    plot(feature_points_c1_mirror{i}(:,1),feature_points_c1_mirror{i}(:,2),'.r')
    pause(0.5);
end
close all

%% Combinación de mascarás (Rectangulo y muestreo de C2)
close all
mask_ROI2={};
mask_final2=mask2;
feature_points_c2_mirror={};
figure;
for i=1:length(mask_final)
    mask_ROI2{i}=mask_final2{i}.*mask_points2;
    [ycoord, xcoord] = find(mask_ROI2{i});
    feature_points_c2_mirror{i}=[xcoord, ycoord];
    subplot(121)
    imshow(texture_img2{i})
    hold on
    plot(feature_points_c2_mirror{i}(:,1),feature_points_c2_mirror{i}(:,2),'.r')
    subplot(122)
    imshow(mask_final2{i})
    hold on
    plot(feature_points_c2_mirror{i}(:,1),feature_points_c2_mirror{i}(:,2),'.r')
    pause(0.5);
end
close all