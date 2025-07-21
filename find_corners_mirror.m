% mask is a logical mask of a squared mirror
% width is the image width.
% mirror_corners are the four corners of the mirror
% By: Eberto Benjumea
% e-mail: ebenjumea@utb.edu.co


function mirror_corners = find_corners_mirror(mask,width)
    % Buscando esquinas superior izq e inferior der.
    corners= detectHarrisFeatures(mask,'MinQuality', 0.01,'FilterSize', 5);
    corndata=corners.selectStrongest(100);
    corn_x=corndata.Location(:,1);
    corn_y=corndata.Location(:,2);
    distance=[];
    for k=1:size(corn_x,1)
        distance(k)=norm([corn_x(k) corn_y(k)]);
    end
    [M,idx_max] = max(distance);
    [M,idx_min] = min(distance);
    corner_sup_izq=[corn_x(idx_min); corn_y(idx_min)];
    corner_inf_der=[corn_x(idx_max); corn_y(idx_max)];
    cor1=[corner_sup_izq corner_inf_der];
    %return mirror_corners

    %Buscando esquinas inferior izq y superior der.
    H=[1 0 -width; 0 1 0; 0 0 1];
    corn_trasl=H*[corn_x';corn_y'; ones(1,size(corn_x,1))];
    corn_trasl_x=corn_trasl(1,:);
    corn_trasl_y=corn_trasl(2,:);
    distance=[];
    for k=1:size(corn_trasl,2)
        distance(k)=norm([corn_trasl_x(k) corn_trasl_y(k)]);
    end  
    [M,idx_max] = max(distance);
    [M,idx_min] = min(distance);
    corner_sup_der=[corn_trasl_x(idx_min); corn_trasl_y(idx_min)];
    corner_inf_izq=[corn_trasl_x(idx_max); corn_trasl_y(idx_max)];
    cor_new_sist=[corner_sup_der corner_inf_izq];
    %Convertimos al sistema original
    cor2_orig=inv(H)*[cor_new_sist; ones(1,2)];
    mirror_corners=[cor1 cor2_orig(1:2,:)];
    mirror_corners=[mirror_corners(:,1) mirror_corners(:,3) mirror_corners(:,2) mirror_corners(:,4)];
end