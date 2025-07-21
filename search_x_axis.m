%This function searchs two points on x-axis in a points cloud
% x, y an z are vectors of the 3D coordinates in camera system
% By: Eberto Benjumea Mendoza
% e-mail: ebenjumea@utb.edu.co
function [ps]=search_x_axis(x,y,z)
% shape of ps=[p0 p1]
ps=[];
for n=1:2
    if n==1
        signo=1;
    else
        signo=-1;
    end
    H=[1 0 -200*signo; 0 1 -200; 0 0 1];
    coor_trasl=H*[x';y'; ones(1,size(x,1))];
    coor_trasl_x=coor_trasl(1,:);
    coor_trasl_y=coor_trasl(2,:);
    distance=[];
    for k=1:size(coor_trasl,2)
        distance(k)=norm([coor_trasl_x(k) coor_trasl_y(k)]);
    end  
    %[M,idx_max] = max(distance);
    [M,idx_min] = min(distance);
    ps(:,n)=[x(idx_min); y(idx_min); z(idx_min)];

end
end