function [mask,final_corners]=find_mirror_nmask(img1, n_corners,v_result)
    
% corner_contour_store=[];
% corner_store=[];
% c_semifinal_corners=[];
% final_corners=[];
img= imadjust(img1);
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
if v_result==1
    subplot(121)
    imshow(img1)
    hold on
    plot(final_corners(:,1),final_corners(:,2),'+b');
    subplot(122) 
    imshow(mask,[])
end



