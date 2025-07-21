function obtain_texture_interpol6_withdatapython(Direc, file,Kt, kt_dist,K,K_dist,Rt,Tt, XcM,YcM,ZcM,img_ir,texture_img,i,height,width,method)
h_ir=size(img_ir,1);
w_ir=size(img_ir,2);




% Kt=params.CameraParameters2.IntrinsicMatrix';
% radial=params.CameraParameters2.RadialDistortion;
% if length(radial)==3
%     k_dist=[radial(1) radial(2) params.CameraParameters2.TangentialDistortion radial(3) ];
% elseif length(radial)==2
%     k_dist=[radial(1) radial(2) params.CameraParameters2.TangentialDistortion  0];
% end
% 
% R=params.RotationOfCamera2';
% T=params.TranslationOfCamera2';


% ojo
% R=params.RotationOfCamera2';
% eul = rotm2eul(R);
% ang=eul*180/pi;
% ang(2)=-2.2;
% ang=ang*pi/180;
% R = eul2rotm(ang);

%%%%%

M_c1_to_t=[Rt Tt; 0 0 0 1];

x=XcM';
y=YcM';
z=ZcM';

% disp(['x has ' num2str(sum(sum(isnan(x)))) ' NaN data'])
% disp(['y has ' num2str(sum(sum(isnan(y)))) ' NaN data'])
% disp(['z has ' num2str(sum(sum(isnan(z)))) ' NaN data'])
%---------Convert points from VIS coordinates to IR camera coordinates
p=[x; y; z; ones(size(x))];
pt=M_c1_to_t*p;
pt=pt(1:3,:);
% Normalizamos los puntos (Normalized image coordinate) 
% pt=pt./pt(3,:);
% coord=pt(1:2,:);



% ---------------------------------------------

% ----Elije entre estas dos----Normalizamos los puntos (Normalized image coordinate) 
% pt=pt./pt(3,:);
% coord=pt(1:2,:);
%----------------------------
pixel_t=Kt*pt;
pixel_t=pixel_t./pixel_t(3,:);
norm_pixel_t=pinv(Kt)*pixel_t;
coord=norm_pixel_t(1:2,:);





%-------------------------------
%- --------Add distortions and convert to image coordinates------------------
%pixel_distun=distortion_points(Kt,k_dist,coord);

%pixel_distun = comp_distortion_oulu(coord,k_dist);
%pixel_distun=Kt*[pixel_distun; ones(1,size(pixel_distun,2))];
pixel_distun=distortion_points(Kt,kt_dist,coord);
%------------------------------------------------------------------------

%pixel_t=Kt*pt;
pixel_distun=pixel_distun./pixel_distun(3,:);
%pixel_distun=pt./pt(3,:);
%pixel_distun=Kt*pixel_distun;
% figure;
% subplot(121)
% imshow(zeros(300,300))
% hold on
% plot(pixel_t(1,:),pixel_t(2,:),'.r')
% line([0 256],[192 192]);
% line([256 256],[0 192]);
% title('Without distortion')
% subplot(122)
% imshow(zeros(300,300))
% hold on
% plot(pixel_distun(1,:),pixel_distun(2,:),'.r')
% line([0 256],[192 192]);
% line([256 256],[0 192]);
% title('With distortion')
% x=XcM';
% y=YcM';
% z=ZcM';

% disp(['x has ' num2str(sum(sum(isnan(x)))) ' NaN data'])
% disp(['y has ' num2str(sum(sum(isnan(y)))) ' NaN data'])
% disp(['z has ' num2str(sum(sum(isnan(z)))) ' NaN data'])
%---------Convert points from VIS coordinates to IR camera coordinates
%pt=[x'; y'; z'];
% pt=M_c1_to_t*p;
% pt=p(1:3,:);
% Normalizamos los puntos (Normalized image coordinate) 
% pt=pt./pt(3,:);
% coord=pt(1:2,:);

% aa=50;
% figure(900),
% plot3(pt(1,1:aa:end),pt(2,1:aa:end),pt(3,1:aa:end),'.')
% graph_cam(M_c1_to_t,1)
% graph_cam([1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1],2)
% axis equal
% 
% figure(901)
% plot3(x(1:aa:end),y(1:aa:end),z(1:aa:end),'.')
% graph_cam([1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1],1)
% graph_cam([R t],2)
% axis equal


% ---------------------------------------------

% ----Elije entre estas dos----Normalizamos los puntos (Normalized image coordinate) 
% pt=pt./pt(3,:);
% coord=pt(1:2,:);
%----------------------------
% pixel_t=Kt*pt;
% pixel_t=pixel_t./pixel_t(3,:);
% norm_pixel_t=pinv(Kt)*pixel_t;
% coord=norm_pixel_t(1:2,:);


%imagePoints = worldToImage(params.CameraParameters2,eye(3),zeros(1,3),pt');


%-------------------------------
%- --------Add distortions and convert to image coordinates------------------
%pixel_distun=distortion_points(Kt,k_dist,coord);

%pixel_distun = comp_distortion_oulu(coord,k_dist);
%pixel_distun=Kt*[pixel_distun; ones(1,size(pixel_distun,2))];
%pixel_distun=distortion_points(Kt,kt_dist,coord);
%------------------------------------------------------------------------
%pixel_distun=[coord; ones(1,size(coord,2))];
%pixel_t=Kt*pt;
%pixel_distun=pixel_distun./pixel_distun(3,:);
%pixel_distun=pt./pt(3,:);
%pixel_distun=Kt*pixel_distun;
% figure;
% subplot(121)
% imshow(zeros(300,300))
% hold on
% plot(pixel_t(1,:),pixel_t(2,:),'.r')
% line([0 256],[192 192]);
% line([256 256],[0 192]);
% title('Without distortion')
% subplot(122)
% imshow(zeros(300,300))
% hold on
% plot(pixel_distun(1,:),pixel_distun(2,:),'.r')
% line([0 256],[192 192]);
% line([256 256],[0 192]);
% title('With distortion')

% ----------------Interpolamos para buscar textura a nivel subpixel--------
[x_ir, y_ir]=meshgrid(1:size(img_ir,2),1:size(img_ir,1));

% x_ir=x_ir(:)';
% y_ir=y_ir(:)';
% aux_img_ir=img_ir(:)';
% mask_pixel_distun=isnan(pixel_distun);
% pixel_distun(mask_pixel_distun)=0;
% texture_ir= griddata(x_ir,y_ir,aux_img_ir,pixel_distun(1,:),pixel_distun(2,:),"cubic");


% x_ir=x_ir(:)';
% y_ir=y_ir(:)';
% aux_img_ir=img_ir(:)';
aux_img_ir=img_ir;
mask_pixel_distun=~isnan(pixel_distun(1,:));
% pixel_distun(isnan(pixel_distun))=[];
%texture_ir= griddata(x_ir,y_ir,aux_img_ir,pixel_distun(1,mask_pixel_distun),pixel_distun(2,mask_pixel_distun),"cubic");
texture_ir= interp2(x_ir,y_ir,aux_img_ir,pixel_distun(1,:),pixel_distun(2,:),"cubic");
figure;
subplot(131)
plot(pixel_distun(1,:))
subplot(132)
plot(pixel_distun(2,:))
subplot(133)
plot(texture_ir)

figure, imagesc(img_ir)
hold on
plot(pixel_distun(1,mask_pixel_distun),pixel_distun(2,mask_pixel_distun),'.r')


% pixel_distun_round=floor([pixel_distun(1,:);pixel_distun(2,:)]);
% figure, imagesc(img_ir)
% hold on
% plot(pixel_distun_round(1,:),pixel_distun_round(2,:),'.r')


%%%%%%%%%%%%%%% Ojo
pixel_distun_round=[pixel_distun(1,:);pixel_distun(2,:)];
for kk=1:size(pixel_distun_round,2)
    %kk
    if isnan(pixel_distun_round(1,kk)) || isnan(pixel_distun_round(2,kk))
        texture_round(kk)=NaN;
        
    elseif abs(pixel_distun_round(1,kk))>256 || abs(pixel_distun_round(2,kk))>192
        texture_round(kk)=NaN;

    elseif pixel_distun_round(1,kk)<1 || pixel_distun_round(2,kk)<1
        texture_round(kk)=NaN;

    elseif pixel_distun_round(1,kk)<=256 && pixel_distun_round(2,kk)<=192
        texture_round(kk)=texture_ir(kk);
    else
        pixel_distun_round(2,kk)
        pixel_distun_round(1,kk)
    end

end

%     [X,Y] = meshgrid(1:size(Itex,2), 1:size(Itex,1));
%     PhasexPoints = interp2(X,Y,fase_c1{i,1},puntos(:,1),puntos(:,2));

%%%%%%%%%%%%%%
% for kk=1:size(pixel_distun_round,2)
%     %kk
%     if isnan(pixel_distun_round(1,kk)) || isnan(pixel_distun_round(2,kk))
%         texture_round(kk)=NaN;
%         
%     elseif abs(pixel_distun_round(1,kk))>256 || abs(pixel_distun_round(2,kk))>192
%         texture_round(kk)=NaN;
% 
%     elseif pixel_distun_round(1,kk)<1 || pixel_distun_round(2,kk)<1
%         texture_round(kk)=NaN;
% 
%     elseif pixel_distun_round(1,kk)<=256 && pixel_distun_round(2,kk)<=192
%         texture_round(kk)=img_ir(pixel_distun_round(2,kk),pixel_distun_round(1,kk));
%     else
%         pixel_distun_round(2,kk)
%         pixel_distun_round(1,kk)
%     end
% end
%a=1
% parpool('local')
% options = optimoptions('solvername','UseParallel',true);
% delete(gcp)
[x_hat, y_hat,x_t,y_t,z_t, ir_texture,nan_idx]=textura_ir_cercano_parfor2(pixel_distun(1,:),pixel_distun(2,:),pt(1,:),pt(2,:),pt(3,:),texture_round,img_ir);
x_hat=pixel_distun(1,:);
y_hat=pixel_distun(2,:);
x_t=pt(1,:);
y_t=pt(2,:);
z_t=pt(3,:);
ir_texture=texture_round;
% x_hat=x_hat(~nan_idx);
% y_hat=y_hat(~nan_idx);
% x_t=x_t(~nan_idx);
% y_t=y_t(~nan_idx);
z_t(nan_idx)=NaN;
ir_texture(nan_idx)=NaN;


figure, imagesc(img_ir)
hold on
plot(x_hat,y_hat,'.r')

%scatter3(pixel_distun(1,mask_pixel_distun),pixel_distun(2,mask_pixel_distun),pixel_distun(2,mask_pixel_distun)*0,1,texture_ir);

M_t_to_c1=pinv(M_c1_to_t);
p_t_filtered=[x_t; y_t; z_t; ones(1,length(x_t))];

K1=K;
k1_dist=K_dist;
% K1=params.CameraParameters1.IntrinsicMatrix';
% radial=params.CameraParameters1.RadialDistortion;
% if length(radial)==3
%     k1_dist=[radial(1) radial(2) params.CameraParameters1.TangentialDistortion radial(3) ];
% elseif length(radial)==2
%     k1_dist=[radial(1) radial(2) params.CameraParameters1.TangentialDistortion];
% end

p_c1_filtered=M_t_to_c1*p_t_filtered;
pixel_c1_filtered=K1*p_c1_filtered(1:3,:);
pixel_c1_filtered=pixel_c1_filtered./pixel_c1_filtered(3,:);
% Buscamos la ubicacion real de los pixeles teniendo en cuenta las
% distorsiones
norm_pixel_c1=pinv(K1)*pixel_c1_filtered;
coord_c1=norm_pixel_c1(1:2,:);
pixel_dist=distortion_points(K1,k1_dist,coord_c1);
pixel_dist=pixel_dist./pixel_dist(3,:);
pixel_c1_filtered=pixel_dist;

%figure, imshow(zeros(height,width))
figure, imshow(texture_img)
hold on
plot(pixel_c1_filtered(1,:),pixel_c1_filtered(2,:),'.r')
%texture_ir_final=ir_texture;
% 
p_c1=[x; y; z];
pixel_c1=K1*p_c1;
pixel_c1=pixel_c1./pixel_c1(3,:);
% Buscamos la ubicacion real de los pixeles teniendo en cuenta las
% distorsiones
norm_pixel_c1=pinv(K1)*pixel_c1;
coord_c1=norm_pixel_c1(1:2,:);
pixel_dist=distortion_points(K1,k1_dist,coord_c1);
pixel_dist=pixel_dist./pixel_dist(3,:);
pixel_c1=pixel_dist;


mask_pixel=~isnan(pixel_c1(1,:));
%texture_ir_final= griddata(pixel_c1_filtered(1,:),pixel_c1_filtered(2,:),ir_texture,pixel_c1(1,mask_pixel),pixel_c1(2,mask_pixel),"cubic");
%texture_ir_final = interp2(pixel_c1_filtered(1,:),pixel_c1_filtered(2,:),ir_texture,pixel_c1(1,mask_pixel),pixel_c1(2,mask_pixel),'cubic');

%ir_texture
mask_filtered=~isnan(pixel_c1_filtered(1,:));
% 
phi1x=scatteredInterpolant(pixel_c1_filtered(1,mask_filtered)',pixel_c1_filtered(2,mask_filtered)',ir_texture(mask_filtered)','natural','nearest');
texture_ir_final=phi1x(pixel_c1(1,mask_pixel)',pixel_c1(2,mask_pixel)');
% phi1x=scatteredInterpolant(pixel_c1_filtered',pixel_c1_filtered',ir_texture','natural','nearest');
% texture_ir_final=phi1x(pixel_c1(1,:)',pixel_c1(2,:)');


if i<10
    aux=['0' num2str(i)];
else
    aux=num2str(i);
end
%ptCloud = pointCloud([x(mask_pixel_distun) y(mask_pixel_distun) -z(mask_pixel_distun)],Intensity=texture_ir');
switch method
    case 1
        % Two Stages
        folder_m='Two Stages';
        if exist([Direc,file ,'Two Stages'], 'dir') ~= 7
            mkdir([Direc,file ], 'Two Stages')
        end
    case 2
        % Two Stages 2
        folder_m='Two Stages 2';
        if exist([Direc,file ,'Two Stages 2'], 'dir') ~= 7
            mkdir([Direc,file ], 'Two Stages 2')
        end    
    case 3
        % Digital Features 
        folder_m='Digital Features CITEDI';
        if exist([Direc, file,'Digital Features CITEDI'], 'dir') ~= 7
            mkdir([Direc, file], 'Digital Features CITEDI')
        end 

    case 4
        % Digital Features 
        folder_m='Digital Features 1 UTB';
        if exist([Direc, file,'Digital Features 1 UTB'], 'dir') ~= 7
            mkdir([Direc, file], 'Digital Features 1 UTB')
        end 
    case 5
        % Two Stages 3
        folder_m='Digital Features 2 UTB';
        if exist([Direc,file ,'Digital Features 2 UTB'], 'dir') ~= 7
            mkdir([Direc,file ], 'Digital Features 2 UTB')
        end    
end

texture_ir_final_matrix=ones(size(XcM))*NaN;
texture_ir_final_matrix(mask_pixel')=texture_ir_final;
% XcM=reshape(XcM,[1024, 1280]);
% YcM=reshape(YcM,[1024, 1280]);
% ZcM=reshape(ZcM,[1024, 1280]);
% texture_ir_final_matrix=reshape(texture_ir_final_matrix,[1024, 1280]);
% nube=zeros(1024,1280,3);
% nube(:,:,1)=XcM;
% nube(:,:,2)=YcM;
% nube(:,:,3)=ZcM;
%ptCloud = pointCloud([XcM(mask_pixel') YcM(mask_pixel') -ZcM(mask_pixel')],Intensity=texture_ir_final);
%ptCloud = pointCloud([XcM(mask_pixel') YcM(mask_pixel') ZcM(mask_pixel')],Intensity=texture_ir_final);
%ptCloud = pointCloud([XcM YcM ZcM],Intensity=texture_ir_final_matrix);
ptCloud = pointCloud([XcM YcM ZcM],Intensity=texture_ir_final_matrix);

%ptCloud = pointCloud(nube,Intensity=texture_ir_final_matrix);


pcwrite(ptCloud,[Direc file folder_m '\' 'new_obj_vis_ir_' aux '.ply'],PLYFormat="binary");
pc = pcread([Direc file folder_m '\' 'new_obj_vis_ir_' aux '.ply']);
figure;
pcshow(pc,'BackgroundColor','w');
xlabel('x (mm)')
ylabel('y (mm)')
zlabel('z (mm)')
c = colorbar; 
c.Label.String = '° C';
            view(45,45)
            view(165,60)
%            zt = get(gca, 'ZTick');
%            set(gca, 'ZTick',zt, 'ZTickLabel',fliplr(zt))
set(gca, 'Zdir', 'reverse')
set(gca, 'Xdir', 'reverse')
% x_c=XcM(mask_pixel');
% y_c=YcM(mask_pixel');
% z_c=ZcM(mask_pixel');

x_c=XcM;
y_c=YcM;
z_c=ZcM;
mask_final=mask_pixel';
save([Direc file folder_m '\' 'obj_vis_ir_' aux '.mat'], 'x_c','y_c','z_c','texture_ir_final','mask_final')


%ptCloud = pointCloud([x_t' y_t' -z_t'],Intensity=ir_texture');
ptCloud = pointCloud([x_t' y_t' z_t'],Intensity=ir_texture');
pcwrite(ptCloud,[Direc file folder_m '\' 'new_obj_vis_ir_t_' aux '.ply'],PLYFormat="binary");
pc = pcread([Direc file folder_m '\' 'new_obj_vis_ir_t_' aux '.ply']);
figure;

pcshow(pc,'BackgroundColor','w');
xlabel('x (mm)')
ylabel('y (mm)')
zlabel('z (mm)')
c = colorbar; 
c.Label.String = '° C';

x_t=x_t';
y_t=y_t';
z_t=z_t';
save([Direc file folder_m '\' 'new_obj_vis_ir_t_' aux '.mat'], 'x_t','y_t','z_t','texture_ir_final')


end