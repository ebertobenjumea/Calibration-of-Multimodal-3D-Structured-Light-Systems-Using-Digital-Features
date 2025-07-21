function pixel_distun=distortion_points(k,k_dist,points)
% k: camera matrix
% k_dist: distortions parameters [k1 k2 p1 p2 k3] ó [k1 k2 p1 p2 k3]
% points: normalized points [x/z; y/z] 
coord=points;
r2=coord(1,:).^2+coord(2,:).^2;
r4=r2.^2;
r6=r2.^3;
xy=coord(1,:).*coord(2,:);
if length(k_dist)==5
    pixel_dist=((1+k_dist(1).*r2+k_dist(2).*r4+k_dist(5).*r6).*coord)+[
    (2*k_dist(3).*xy)+(k_dist(4).*(r2+2*(coord(1,:)).^2));
    k_dist(3).*(r2+2*(coord(2,:)).^2)+2*k_dist(4).*xy
    ];
elseif length(k_dist)==4
    pixel_dist=(1+k_dist(1).*r2+k_dist(2).*r4).*coord+[
    2*k_dist(3).*xy+k_dist(4).*(r2+2*(coord(1,:)).^2);
    k_dist(3).*(r2+2*(coord(2,:)).^2)+2*k_dist(4).*xy
    ];
end
%pixel_dist=coord;
% -----------Convertimos a coordenadas imagen-------
pixel_distun=k*[pixel_dist; ones(1,size(pixel_dist,2))];



end