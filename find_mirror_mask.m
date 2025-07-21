function [mask,mask_final,corner_mirror]=find_mirror_mask(img, delta,width,height,group_out,v_result) 

level = graythresh(img);        %Umbralizamos
mask1 = imbinarize(img,level);   
[mask_w_w,numWhite] = bwlabel(mask1);
mask_w_w = bwareaopen(mask_w_w,group_out);  %Retirando manchas blancas
[L,numBlack] = bwlabel(~mask_w_w); %Etiquetamos los agrupamientos de pixeles.
mask= ~bwareaopen(L,group_out);   %Retirando manchas negras
corner_mirror=find_corners_mirror(mask,width);     %Detectando esquinas del espejo
% OJO DISMINUYENDO AREA DE ESPEJO
%delta=0;
corner_mirror(1,1)=corner_mirror(1,1)+delta;
corner_mirror(2,1)=corner_mirror(2,1)+delta;

corner_mirror(1,2)=corner_mirror(1,2)-delta;
corner_mirror(2,2)=corner_mirror(2,2)+delta;

corner_mirror(1,3)=corner_mirror(1,3)-delta;
corner_mirror(2,3)=corner_mirror(2,3)-delta;

corner_mirror(1,4)=corner_mirror(1,4)+delta;
corner_mirror(2,4)=corner_mirror(2,4)-delta;

% OJO DISMINUYENDO AREA DE ESPEJO
mask_final= poly2mask(double(corner_mirror(1,:)),double(corner_mirror(2,:)),height,width);
    
if v_result==1
    subplot(131)
    imshow(img)
    subplot(132) 
    imshow(mask,[])
    hold on
    plot(corner_mirror(1,:),corner_mirror(2,:),'+g')
    plot([corner_mirror(1,:) corner_mirror(1,1)],[corner_mirror(2,:) corner_mirror(2,1)],'r')
    subplot(133)
    imshow(mask_final)
end


end