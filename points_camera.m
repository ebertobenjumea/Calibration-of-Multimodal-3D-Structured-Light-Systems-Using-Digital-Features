
function [screen_points,mask_points]=points_camera(aux_img,grid_size)


mask_points=zeros(size(aux_img));
cnt_points=0;
screen_points=[];
for x=1:size(aux_img,2)
    residuo_x = mod(x,grid_size);
    for y=1:size(aux_img,1)
        residuo_y=mod(y,grid_size);
        if residuo_x==0 && residuo_y==0
            mask_points(y,x)=1;
            cnt_points=cnt_points+1;
            screen_points(cnt_points,:)=[x,y];
        end
    end
end
mask_points=boolean(mask_points);
end