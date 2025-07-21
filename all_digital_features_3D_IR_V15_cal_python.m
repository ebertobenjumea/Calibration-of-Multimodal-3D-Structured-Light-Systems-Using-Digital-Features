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

%% Obtención de fase
disp('Phase obtaining.')
disp(' ')
mode_phase=input(['If you want to obtain the phase from the images, please press 1\n' ...
    'If you want to read the phase from previously processed .mat files, \nPress 2\n Your answer:']);
if mode_phase==1
    disp('Obtaining phase from images...')
    NStep = 18;
    NBits = 7;
    inicio=1;
    carp=camera1;
    fase_c1=estimate_phases(direc,carp,NStep,NBits,n_poses,width,height,format,inicio);
    carp=camera2;
    fase_c2=estimate_phases(direc,carp,NStep,NBits,n_poses,width,height,format,inicio);
    disp('Finished!')
elseif mode_phase==2
    disp('Reading phase from .mat files...')
    fase_c1={};
    fase_c2={};
    for i=1:n_poses
        if i<11, aux=['0' num2str(i-1)]; else, aux=num2str(i-1); end
        load([direc camera1 '\camera 1_pose_' aux '_phases.mat'])
        fase_c1{i,1}=Fx;
        fase_c1{i,2}=Fy;
        load([direc camera2 '\camera 2_pose_' aux '_phases.mat'])
        fase_c2{i,1}=Fx;
        fase_c2{i,2}=Fy;
    end
    disp('Finished!')    
end

%Size of fase_c2: (n_poses, fx fy)
%         fx      fy   
%         .       .
%         .       .   
%         .       .
%     n_poses   n_poses
figure;
for i=1:n_poses
%     fase_c1{i,1}(~mask_final{i})=NaN;
%     fase_c1{i,2}(~mask_final{i})=NaN;
%     fase_c2{i,1}(~mask_final2{i})=NaN;
%     fase_c2{i,2}(~mask_final2{i})=NaN;
    title(['Pose ' num2str(i-1)])
    subplot(221)
    imagesc(fase_c1{i,1})
    axis equal
    subplot(222)
    imagesc(fase_c1{i,2})
    axis equal
    subplot(223)
    imagesc(fase_c2{i,1})
    axis equal
    subplot(224)
    imagesc(fase_c2{i,2})
    axis equal
    sgtitle(['Pose ' num2str(i-1)])
    pause(1);
end

%% Obtenemos fase en puntos de C1 con resolución subpixel 
%    y buscamos correspondencia en proyector vía fase
proyector_points={};
if inst_ask==1
    %W=912; H=1140;
        W=1280; H=800;
elseif inst_ask==2
    W=1280; H=800;
end
P=18; 
figure;
f1_sampled={};

for i=1:n_poses
    tic
    disp(['Pose ' num2str(i-1)])
    % Load point coordinates
    x=feature_points_c1_mirror{i}(:,1);
    y=feature_points_c1_mirror{i}(:,2);
    % Load phase in camera 1
    f1x=fase_c1{i,1};
    f1y=fase_c1{i,2};
    % Search the phase in the points
    phi1x=[];
    phi1y=[];
    for m=1:length(x)
        phi1x(m,1)=f1x(y(m),x(m));
        phi1y(m,1)=f1y(y(m),x(m));
    end

%     f1x=fase_c1{i,1};
%     f1y=fase_c1{i,2};
%     [x,y]=meshgrid(1:width,1:height);
%     phi1x=f1x(:);
%     phi1y=f1y(:);
%     x=x(:);
%     y=y(:);

    coord_x=feature_points_c1_mirror{i}(:,1);
    coord_y=feature_points_c1_mirror{i}(:,2);
    %Interpolamos fase en x
    f1x_int= griddata(double(x),double(y),phi1x,double(coord_x),double(coord_y),"cubic");

%     phi1x=scatteredInterpolant(double(x),double(y),phi1x,'natural','nearest');
%     f1x_int=phi1x(double(coord_x),double(coord_y));
    %Interpolamos fase en y
    f1y_int= griddata(double(x),double(y),phi1y,double(coord_x),double(coord_y),"cubic");
    
%     phi1y=scatteredInterpolant(double(x),double(y),phi1y,'natural','nearest');
%     f1y_int=phi1y(double(coord_x),double(coord_y));
    % figure;
    % plot(f1y_sampled,'.b')
    % hold on
    % plot(f1y_int,'.r')
    f1_sampled{i,1}=f1x_int;
    f1_sampled{i,2}=f1y_int;
    proyector_points{i}(:,1)=P*f1x_int/(2*pi);
    proyector_points{i}(:,2)=P*f1y_int/(2*pi);

    imshow(zeros(H,W))
    hold on
    plot(proyector_points{i}(:,1),proyector_points{i}(:,2),'.w')
    toc
    pause(1);
end

%% Interpolacion griddata

close all
clc

figure;
coor_c2={};
for i=1:n_poses
    % Load points on C1 and C2 sensors
    x=feature_points_c1_mirror{i}(:,1);
    y=feature_points_c1_mirror{i}(:,2);




    % Load phases on x and y for each camera
    phi1x=f1_sampled{i,1};
    phi1y=f1_sampled{i,2};

    x2=feature_points_c2_mirror{i}(:,1);
    y2=feature_points_c2_mirror{i}(:,2);
    f2x=fase_c2{i,1};
    f2y=fase_c2{i,2};
    [row, col]=find(mask_final2{i}==1);
    %figure; imshow(mask_final2{i})
    x2=col;
    y2=row;
    phi2x=[];
    phi2y=[];
    for m=1:length(x2)
        phi2x(m)=f2x(y2(m),x2(m));
        phi2y(m)=f2y(y2(m),x2(m));
    end

%     f2x=fase_c2{i,1};
%     f2y=fase_c2{i,2};
%     [x2,y2]=meshgrid(1:width,1:height);
%     phi2x=f2x(:);
%     phi2y=f2y(:);
%     x2=x2(:);
%     y2=y2(:);


    f1= griddata(phi2x',phi2y',double(x2'),phi1x',phi1y',"cubic");
    f2 = griddata(phi2x',phi2y',double(y2'),phi1x',phi1y',"cubic");
    
    %[Xq,Yq,vq] =griddata(phi2x,phi2y,double(vc2'),phi1x',phi1y',"cubic");
    subplot(121)
    imagesc(fase_c1{i,2})
    %imshow(mask_final{i})
    axis equal
    hold on
    plot(x,y,'.r')
    ylim([0 height])
    hold off
    title(['Pose ' num2str(i-1) ' - Camera 1'])
    subplot(122)
    imagesc(fase_c2{i,2})
    %imshow(mask_final2{i})
    axis equal
    hold on
    plot(f1,f2,'.r')
    ylim([0 height])
    title(['Pose ' num2str(i-1) ' - Camera 2'])
    coor_c2{i}=[f1' f2'];
    pause(1);

end



%% Evaluando precisión
close all
if inst_ask==1
    load('D:\R3D\Adquisiciones\Experimento 54\Digital features\Stereo calibration\calib_stereo_python.mat')
    %load('D:\R3D\Adquisiciones\Experimento 54\Digital features\Stereo calibration\Calib_stereoparams0_081120242.mat')
elseif inst_ask==2
    load('D:\R3D\2023\multimodal\Experimento 44\Digital Features 2\Stereo Calibration\Calib_stereoparams0_270620233.mat')
    %load('D:\R3D\2023\multimodal\Experimento 18\calibracion stereo\Calib_stereoparams0_270620233.mat')
end
K1
K1_dist
K2
K2_dist
R=R_st
t=t_st
E 
F






% st=StereoParams;
% K1=st.CameraParameters1.IntrinsicMatrix'
% dist1=st.CameraParameters1.RadialDistortion
% dist2=st.CameraParameters1.TangentialDistortion
% K2=st.CameraParameters2.IntrinsicMatrix';
% dist21=st.CameraParameters2.RadialDistortion;
% dist22=st.CameraParameters2.TangentialDistortion
% R=st.RotationOfCamera2;
% t=st.TranslationOfCamera2;
% F=st.FundamentalMatrix;
% E=st.EssentialMatrix;
% %save('calib.mat','t',"R","K1","K2","dist1","dist2","F","E")
% dist_k1=[dist1 0 dist2];


figure();
for nm=1:n_poses
    e=[];
    for np=1:size(coor_c2{nm},1)

        p1=[feature_points_c1_mirror{nm}(np,:) 1];
        p2=[coor_c2{nm}(np,:) 1];
        e(np)=abs(p1*F*p2');
    end   
    plot(e)
    title(['Error absoluto para p1.T*F*p2 para pose ' num2str(nm-1)])
    pause(1);
end

%%
close all
clc
%st=StereoParams;
%cam1_params=st.CameraParameters1;
%cam2_params=st.CameraParameters2;
addpath('C:\Users\Eberto Benjumea\Desktop\OneDrive - Universidad Tecnológica de Bolívar\Doctorado\Proyectos\Codigos Reconstruccion 3D Learning\codigos calibración\code\Funciones_Hibrido')
%n_poses=3;
worldPoints={};
model_points={};
planes=[];
n=[];
p0_store=[];
H_stores={};

figure(1);
for i=1:n_poses
    %i=16
    disp(['Pose ' num2str(i-1)])
    mask_c2=~isnan(coor_c2{i});
    num_non_processed=sum(isnan(coor_c2{i}));
    disp(['Number of non processed points: ' num2str(num_non_processed(1))])
    coor_c2_aux=coor_c2{i}(mask_c2(:,1),:);
    coor_c1_aux=feature_points_c1_mirror{i}(mask_c2(:,1),:);
    norm_coor_c1_aux=pinv(K1)*[double(coor_c1_aux'); ones(1,size(coor_c1_aux,1))];
    norm_coor_c2_aux=pinv(K2)*[double(coor_c2_aux'); ones(1,size(coor_c2_aux,1))];
    norm_coor_c1_aux=norm_coor_c1_aux(1:2,:);
    norm_coor_c2_aux=norm_coor_c2_aux(1:2,:);
    coor_c1_undist =comp_distortion_oulu(norm_coor_c1_aux,K1_dist);
    coor_c2_undist =comp_distortion_oulu(norm_coor_c2_aux,K2_dist);
    coor_c1_undist=K1*[coor_c1_undist; ones(1,size(coor_c1_undist,2))];
    coor_c2_undist=K2*[coor_c2_undist; ones(1,size(coor_c2_undist,2))];
    coor_c1_undist=coor_c1_undist(1:2,:);
    coor_c2_undist=coor_c2_undist(1:2,:);
    data=pyrunfile("triangulate.py","X",K1=py.numpy.array(K1), K2=py.numpy.array(K2), R=py.numpy.array(R_st), t=py.numpy.array(t_st),pk1=py.numpy.array(coor_c1_undist),pk2=py.numpy.array(coor_c2_undist));
    %res = pyrunfile("addac.py","z",x=3,y=2)
    data= double(data);
    %coor_c1_undist = undistortPoints(double(coor_c1_aux),cam1_params);
    %coor_c2_undist = undistortPoints(coor_c2_aux,cam2_params);
    %size(coor_c1_undist)
    %size(coor_c2_undist)
    %worldPoints{i} = triangulate(coor_c1_undist,coor_c2_undist,StereoParams);
%     XcM=worldPoints{i}(:,1);
%     YcM=worldPoints{i}(:,2);
%     ZcM=worldPoints{i}(:,3);
    XcM=data(1,:)';
    YcM=data(2,:)';
    ZcM=data(3,:)';
    

%    figure;
    subplot(231)
    plot3(XcM,YcM,ZcM,".")
    view([45 45])
    xlabel('x (mm)')
    ylabel('y (mm)')
    zlabel('z (mm)')
    title('3D reconstructed points using stereo calibration data')
    %axis equal

%     prom_z=mean(z);
%     s_z = std(z);
%     z=z(z<prom_z+2*s_z & z>prom_z-2*s_z);
%     x=x(z<prom_z+2*s_z & z>prom_z-2*s_z);
%     y=y(z<prom_z+2*s_z & z>prom_z-2*s_z);


    x=XcM;
    y=YcM;
    z=ZcM;
    DM = [x, y, ones(size(z))];                             % Design Matrix
    %B = DM\z;                                               % Estimate Parameters
    B=lsqminnorm(DM,z);
    [X,Y] = meshgrid(linspace(min(x),max(x),50), linspace(min(y),max(y),50));
    Z = B(1)*X + B(2)*Y + B(3)*ones(size(X));
    planes(i,:)=[B(1) B(2) B(3)];
    n(:,i)=[B(1); B(2); -1];

%     subplot(232)
%figure;
%     surf(X, Y, Z)
%     hold on
%     plot3(x,y,z,'.k')
%     view([45 45])
%     hold off
%     xlabel('x (mm)')
%     ylabel('y (mm)')
%     zlabel('z (mm)')
%     title('Ideal 3D plane')
    
    subplot(232)
    Z = B(1)*x + B(2)*y + B(3)*ones(size(x));
    plot3(x,y,Z,'.r')
    view([45 45])
    title('Reconstructed ideal 3D plane')
    xlabel('x (mm)')
    ylabel('y (mm)')
    zlabel('z (mm)')
    
    subplot(233)
    error_z=z-Z;
    plot3(x,y,error_z,'.r')
%     dif_p=[x-x y-y z-Z];
%     error=[];
%     for j=1:size(dif_p,1)
%         error(j)=norm(dif_p(j,:));
%     end
%     plot3(x,y,error,'.r')
    view([45 45])
    %zlim([-2 2])
    xlabel('x (mm)')
    ylabel('y (mm)')
    zlabel('z (mm)')
    title('Plane fitting error')
    

    subplot(234)
    x0=x(1);
    y0=y(1);
    z0=Z(1);
    norm_plan=n(:,i);
    nz=norm_plan./norm(norm_plan);
    %figure; 

    % Busqueda de eje x





    x0_point_c1_aux=pinv(K1)*[double(corner_store{i}(1,:)'); 1];
    x0_point_c2_aux=pinv(K2)*[double(corner_store2{i}(1,:)'); 1];
    x0_point_c1_aux=x0_point_c1_aux(1:2,:);
    x0_point_c2_aux=x0_point_c2_aux(1:2,:);
    x0_point_c1_undist =comp_distortion_oulu(x0_point_c1_aux,K1_dist);
    x0_point_c2_undist =comp_distortion_oulu(x0_point_c2_aux,K2_dist);
    x0_point_c1_undist=K1*[x0_point_c1_undist; ones(1,size(x0_point_c1_undist,2))];
    x0_point_c2_undist=K2*[x0_point_c2_undist; ones(1,size(x0_point_c2_undist,2))];
    x0_point_c1_undist=x0_point_c1_undist(1:2,:);
    x0_point_c2_undist=x0_point_c2_undist(1:2,:);
    data=pyrunfile("triangulate.py","X",K1=py.numpy.array(K1), K2=py.numpy.array(K2), R=py.numpy.array(R_st), t=py.numpy.array(t_st),pk1=py.numpy.array(x0_point_c1_undist),pk2=py.numpy.array(x0_point_c2_undist));
    x0_point= double(data)';

    x1_point_c1_aux=pinv(K1)*[double(corner_store{i}(2,:)'); 1];
    x1_point_c2_aux=pinv(K2)*[double(corner_store2{i}(2,:)'); 1];
    x1_point_c1_aux=x1_point_c1_aux(1:2,:);
    x1_point_c2_aux=x1_point_c2_aux(1:2,:);
    x1_point_c1_undist =comp_distortion_oulu(x1_point_c1_aux,K1_dist);
    x1_point_c2_undist =comp_distortion_oulu(x1_point_c2_aux,K2_dist);
    x1_point_c1_undist=K1*[x1_point_c1_undist; ones(1,size(x1_point_c1_undist,2))];
    x1_point_c2_undist=K2*[x1_point_c2_undist; ones(1,size(x1_point_c2_undist,2))];
    x1_point_c1_undist=x1_point_c1_undist(1:2,:);
    x1_point_c2_undist=x1_point_c2_undist(1:2,:);
    data=pyrunfile("triangulate.py","X",K1=py.numpy.array(K1), K2=py.numpy.array(K2), R=py.numpy.array(R_st), t=py.numpy.array(t_st),pk1=py.numpy.array(x1_point_c1_undist),pk2=py.numpy.array(x1_point_c2_undist));
    x1_point= double(data)';


%     x0_point_c1 = undistortPoints(double(corner_store{i}(1,:)'),cam1_params);
%     x0_point_c2 = undistortPoints(double(corner_store2{i}(1,:)),cam2_params);
%     x0_point = triangulate(x0_point_c1,x0_point_c2,StereoParams);
%     x1_point_c1 = undistortPoints(double(corner_store{i}(2,:)),cam1_params);
%     x1_point_c2 = undistortPoints(double(corner_store2{i}(2,:)),cam2_params);
%     x1_point = triangulate(x1_point_c1,x1_point_c2,StereoParams);
%     hold on
%     plot3(x0_point(1,1),x0_point(1,2),x0_point(1,3),'ok')
%     plot3(x1_point(1,1),x1_point(1,2),x1_point(1,3),'ob')
    x0_point(3) = B(1)*x0_point(1)+ B(2)*x0_point(2) + B(3)*ones(size(x0_point(1)));
    x1_point(3) = B(1)*x1_point(1)+ B(2)*x1_point(2) + B(3)*ones(size(x1_point(1)));
    p_0=x0_point';
    p_1=x1_point';




%     ps=search_x_axis(x,y,Z);
%     p_0=ps(:,2);%+[10;-10];
%     p_1=ps(:,1);
    plot(x,y,'.')
    hold on
    plot(p_0(1),p_0(2),'*r')
    plot(p_1(1),p_1(2),'*k')
    hold off
    Nx=p_1-p_0;
    nx=Nx./norm(Nx);
    nz_x=[0 -nz(3) nz(2); nz(3) 0 -nz(1); -nz(2) nz(1) 0];
    ny=nz_x*nx;
    H_c_m=[nx';ny'; nz'];
    H_stores{i}=H_c_m;
    p0_store(:,i)=p_0;
    model_points{i}=H_c_m*[x'-p_0(1);y'-p_0(2);Z'-p_0(3)];
    xlabel('x (mm)')
    ylabel('y (mm)')
    zlabel('z (mm)')
    title('x and y coordinates')
    
    subplot(235)
    plot3(model_points{i}(1,:),model_points{i}(2,:),model_points{i}(3,:),'.k')
    zlim([-1 1])
    xlabel('x (mm)')
    ylabel('y (mm)')
    zlabel('z (mm)')
    title('Points on model coordinate system')
    sgtitle(['Pose ' num2str(i-1)])
    pause(0.5);
%     close 
end

% figure; imshow(texture_img{i}); hold on
% plot(corner_store{i}(5,1),corner_store{i}(5,2),'+b')
% x0_point_c1 = undistortPoints(double(corner_store{i}(6,:)),cam1_params);
% x0_point_c2 = undistortPoints(double(corner_store2{i}(6,:)),cam2_params);
% x0_point = triangulate(x0_point_c1,x0_point_c2,StereoParams);
% x1_point_c1 = undistortPoints(double(corner_store{i}(5,:)),cam1_params);
% x1_point_c2 = undistortPoints(double(corner_store2{i}(5,:)),cam2_params);
% x1_point = triangulate(x1_point_c1,x1_point_c2,StereoParams);
% hold on
% plot3(x0_point(1,1),x0_point(1,2),x0_point(1,3),'ok')
% plot3(x1_point(1,1),x1_point(1,2),x1_point(1,3),'ob')
% 
% 
% x0_point(3) = B(1)*x0_point(1)+ B(2)*x0_point(2) + B(3)*ones(size(x0_point(1)));
% x1_point(3) = B(1)*x1_point(1)+ B(2)*x1_point(2) + B(3)*ones(size(x1_point(1)));
%%
    i=32
    disp(['Pose ' num2str(i-1)])
    mask_c2=~isnan(coor_c2{i});
    num_non_processed=sum(isnan(coor_c2{i}));
    disp(['Number of non processed points: ' num2str(num_non_processed(1))])
    coor_c2_aux=coor_c2{i}(mask_c2(:,1),:);
    coor_c1_aux=feature_points_c1_mirror{i}(mask_c2(:,1),:);
    norm_coor_c1_aux=pinv(K1)*[double(coor_c1_aux'); ones(1,size(coor_c1_aux,1))];
    norm_coor_c2_aux=pinv(K2)*[double(coor_c2_aux'); ones(1,size(coor_c2_aux,1))];
    norm_coor_c1_aux=norm_coor_c1_aux(1:2,:);
    norm_coor_c2_aux=norm_coor_c2_aux(1:2,:);
    coor_c1_undist =comp_distortion_oulu(norm_coor_c1_aux,K1_dist);
    coor_c2_undist =comp_distortion_oulu(norm_coor_c2_aux,K2_dist);
    coor_c1_undist=K1*[coor_c1_undist; ones(1,size(coor_c1_undist,2))];
    coor_c2_undist=K2*[coor_c2_undist; ones(1,size(coor_c2_undist,2))];
    coor_c1_undist=coor_c1_undist(1:2,:);
    coor_c2_undist=coor_c2_undist(1:2,:);
    data=pyrunfile("triangulate.py","X",K1=py.numpy.array(K1), K2=py.numpy.array(K2), R=py.numpy.array(R_st), t=py.numpy.array(t_st),pk1=py.numpy.array(coor_c1_undist),pk2=py.numpy.array(coor_c2_undist));
    %res = pyrunfile("addac.py","z",x=3,y=2)
    data= double(data);
    %coor_c1_undist = undistortPoints(double(coor_c1_aux),cam1_params);
    %coor_c2_undist = undistortPoints(coor_c2_aux,cam2_params);
    %size(coor_c1_undist)
    %size(coor_c2_undist)
    %worldPoints{i} = triangulate(coor_c1_undist,coor_c2_undist,StereoParams);
%     XcM=worldPoints{i}(:,1);
%     YcM=worldPoints{i}(:,2);
%     ZcM=worldPoints{i}(:,3);
    XcM=data(1,:)';
    YcM=data(2,:)';
    ZcM=data(3,:)';
    

%    figure;
    subplot(221)
    plot3(XcM,YcM,ZcM,".")
    view([45 45])
    xlabel('x (mm)')
    ylabel('y (mm)')
    zlabel('z (mm)')
    %axis equal
    %title('3D reconstructed points using stereo calibration data')
    %axis equal

%     prom_z=mean(z);
%     s_z = std(z);
%     z=z(z<prom_z+2*s_z & z>prom_z-2*s_z);
%     x=x(z<prom_z+2*s_z & z>prom_z-2*s_z);
%     y=y(z<prom_z+2*s_z & z>prom_z-2*s_z);


    x=XcM;
    y=YcM;
    z=ZcM;
    DM = [x, y, ones(size(z))];                             % Design Matrix
    %B = DM\z;                                               % Estimate Parameters
    B=lsqminnorm(DM,z);
    [X,Y] = meshgrid(linspace(min(x),max(x),50), linspace(min(y),max(y),50));
    Z = B(1)*X + B(2)*Y + B(3)*ones(size(X));
    planes(i,:)=[B(1) B(2) B(3)];
    n(:,i)=[B(1); B(2); -1];

%     subplot(232)
%figure;
%     surf(X, Y, Z)
%     hold on
%     plot3(x,y,z,'.k')
%     view([45 45])
%     hold off
%     xlabel('x (mm)')
%     ylabel('y (mm)')
%     zlabel('z (mm)')
%     title('Ideal 3D plane')
    
    subplot(222)
    Z = B(1)*x + B(2)*y + B(3)*ones(size(x));
    plot3(x,y,Z,'.r')
    view([45 45])
    %title('Reconstructed ideal 3D plane')
    xlabel('x (mm)')
    ylabel('y (mm)')
    zlabel('z (mm)')
    %axis equal
    
%    subplot(233)
    error_z=z-Z;
 %   plot3(x,y,error_z,'.r')
%     dif_p=[x-x y-y z-Z];
%     error=[];
%     for j=1:size(dif_p,1)
%         error(j)=norm(dif_p(j,:));
%     end
%     plot3(x,y,error,'.r')
%     view([45 45])
%     %zlim([-2 2])
%     xlabel('x (mm)')
%     ylabel('y (mm)')
%     zlabel('z (mm)')
%     title('Plane fitting error')
    

%    subplot(234)
    x0=x(1);
    y0=y(1);
    z0=Z(1);
    norm_plan=n(:,i);
    nz=norm_plan./norm(norm_plan);
    %figure; 

    % Busqueda de eje x





    x0_point_c1_aux=pinv(K1)*[double(corner_store{i}(1,:)'); 1];
    x0_point_c2_aux=pinv(K2)*[double(corner_store2{i}(1,:)'); 1];
    x0_point_c1_aux=x0_point_c1_aux(1:2,:);
    x0_point_c2_aux=x0_point_c2_aux(1:2,:);
    x0_point_c1_undist =comp_distortion_oulu(x0_point_c1_aux,K1_dist);
    x0_point_c2_undist =comp_distortion_oulu(x0_point_c2_aux,K2_dist);
    x0_point_c1_undist=K1*[x0_point_c1_undist; ones(1,size(x0_point_c1_undist,2))];
    x0_point_c2_undist=K2*[x0_point_c2_undist; ones(1,size(x0_point_c2_undist,2))];
    x0_point_c1_undist=x0_point_c1_undist(1:2,:);
    x0_point_c2_undist=x0_point_c2_undist(1:2,:);
    data=pyrunfile("triangulate.py","X",K1=py.numpy.array(K1), K2=py.numpy.array(K2), R=py.numpy.array(R_st), t=py.numpy.array(t_st),pk1=py.numpy.array(x0_point_c1_undist),pk2=py.numpy.array(x0_point_c2_undist));
    x0_point= double(data)';

    x1_point_c1_aux=pinv(K1)*[double(corner_store{i}(2,:)'); 1];
    x1_point_c2_aux=pinv(K2)*[double(corner_store2{i}(2,:)'); 1];
    x1_point_c1_aux=x1_point_c1_aux(1:2,:);
    x1_point_c2_aux=x1_point_c2_aux(1:2,:);
    x1_point_c1_undist =comp_distortion_oulu(x1_point_c1_aux,K1_dist);
    x1_point_c2_undist =comp_distortion_oulu(x1_point_c2_aux,K2_dist);
    x1_point_c1_undist=K1*[x1_point_c1_undist; ones(1,size(x1_point_c1_undist,2))];
    x1_point_c2_undist=K2*[x1_point_c2_undist; ones(1,size(x1_point_c2_undist,2))];
    x1_point_c1_undist=x1_point_c1_undist(1:2,:);
    x1_point_c2_undist=x1_point_c2_undist(1:2,:);
    data=pyrunfile("triangulate.py","X",K1=py.numpy.array(K1), K2=py.numpy.array(K2), R=py.numpy.array(R_st), t=py.numpy.array(t_st),pk1=py.numpy.array(x1_point_c1_undist),pk2=py.numpy.array(x1_point_c2_undist));
    x1_point= double(data)';


%     x0_point_c1 = undistortPoints(double(corner_store{i}(1,:)'),cam1_params);
%     x0_point_c2 = undistortPoints(double(corner_store2{i}(1,:)),cam2_params);
%     x0_point = triangulate(x0_point_c1,x0_point_c2,StereoParams);
%     x1_point_c1 = undistortPoints(double(corner_store{i}(2,:)),cam1_params);
%     x1_point_c2 = undistortPoints(double(corner_store2{i}(2,:)),cam2_params);
%     x1_point = triangulate(x1_point_c1,x1_point_c2,StereoParams);
%     hold on
%     plot3(x0_point(1,1),x0_point(1,2),x0_point(1,3),'ok')
%     plot3(x1_point(1,1),x1_point(1,2),x1_point(1,3),'ob')
    x0_point(3) = B(1)*x0_point(1)+ B(2)*x0_point(2) + B(3)*ones(size(x0_point(1)));
    x1_point(3) = B(1)*x1_point(1)+ B(2)*x1_point(2) + B(3)*ones(size(x1_point(1)));
    p_0=x0_point';
    p_1=x1_point';




%     ps=search_x_axis(x,y,Z);
%     p_0=ps(:,2);%+[10;-10];
%     p_1=ps(:,1);
%     plot(x,y,'.')
%     hold on
%     plot(p_0(1),p_0(2),'*r')
%     plot(p_1(1),p_1(2),'*k')
%     hold off
    Nx=p_1-p_0;
    nx=Nx./norm(Nx);
    nz_x=[0 -nz(3) nz(2); nz(3) 0 -nz(1); -nz(2) nz(1) 0];
    ny=nz_x*nx;
    H_c_m=[nx';ny'; nz'];
    H_stores{i}=H_c_m;
    p0_store(:,i)=p_0;
    model_points{i}=H_c_m*[x'-p_0(1);y'-p_0(2);Z'-p_0(3)];
%     xlabel('x (mm)')
%     ylabel('y (mm)')
%     zlabel('z (mm)')
%     title('x and y coordinates')
    
    subplot(223)
    plot3(model_points{i}(1,:),model_points{i}(2,:),model_points{i}(3,:),'.g')
    zlim([-1 1])
    xlabel('x (mm)')
    ylabel('y (mm)')
    zlabel('z (mm)')
        view([45 45])
        %axis equal
    %title('Points on model coordinate system')
    %sgtitle(['Pose ' num2str(i-1)])

%% Adaptative grid
close all
pattern_size=[15,15];
i=1;
disp(['Pose ' num2str(i-1)])
% obj_x=x_pattern(:);
% obj_y=y_pattern(:);



corners_c1_aux=pinv(K1)*[double(corner_store{i})'; ones(1,size(corner_store{i}',2))];
corners_c2_aux=pinv(K2)*[double(corner_store2{i})'; ones(1,size(corner_store{i}',2))];
corners_c1_aux=corners_c1_aux(1:2,:);
corners_c2_aux=corners_c2_aux(1:2,:);
corners_c1_undist =comp_distortion_oulu(corners_c1_aux,K1_dist);
corners_c2_undist =comp_distortion_oulu(corners_c2_aux,K2_dist);
corners_c1_undist=K1*[corners_c1_undist; ones(1,size(corners_c1_undist,2))];
corners_c2_undist=K2*[corners_c2_undist; ones(1,size(corners_c2_undist,2))];
corners_c1_undist=corners_c1_undist(1:2,:);
corners_c2_undist=corners_c2_undist(1:2,:);
data=pyrunfile("triangulate.py","X",K1=py.numpy.array(K1), K2=py.numpy.array(K2), R=py.numpy.array(R_st), t=py.numpy.array(t_st),pk1=py.numpy.array(corners_c1_undist),pk2=py.numpy.array(corners_c2_undist));
corners_3d= double(data)';


%corners_c1 = undistortPoints(double(corner_store{i}),cam1_params);
%corners_c2 = undistortPoints(double(corner_store2{i}),cam2_params);
%corners_3d = triangulate(corners_c1,corners_c2,StereoParams);
B=planes(i,:);

C = B(1)*corners_3d(:,1) + B(2)*corners_3d(:,2)  + B(3)*ones(size(corners_3d(:,1)));
corners_3d_ideal=[corners_3d(:,1) corners_3d(:,2) C];
figure;
plot3(corners_3d_ideal(:,1),corners_3d_ideal(:,2),corners_3d_ideal(:,3),'ok')
xlabel('x')
ylabel('y')
zlabel('z')

corners_3d_model=H_stores{i}*(corners_3d_ideal'-corners_3d_ideal(6,:)');

figure;
plot3(corners_3d_model(1,:),corners_3d_model(2,:),corners_3d_model(3,:),'ok')
xlabel('x')
ylabel('y')
zlabel('z')
zlim([-1 1])

%y_axis=[corners_3d_model(:,6) corners_3d_model(:,1)];
%x_axis=[corners_3d_model(:,6) corners_3d_model(:,5)];
x_distance=norm([corners_3d_model(:,1)-corners_3d_model(:,2)]);
y_distance=norm([corners_3d_model(:,1)-corners_3d_model(:,6)]);
b61=corners_3d_model(:,6)-corners_3d_model(:,1);
a31=corners_3d_model(:,3)-corners_3d_model(:,1);
y_distance_plus=norm(a31-(dot(a31,b61)/dot(b61,b61)*b61));



% ps=search_x_axis(model_points{i}(1,:)',model_points{i}(2,:)',model_points{i}(3,:)');
% pd=search_y_axis(model_points{i}(1,:)',model_points{i}(2,:)',model_points{i}(3,:)');
% corners_3d=[ps pd];
% for j=1:4
%     distance(j)=norm([corners_3d(:,1)-corners_3d(:,j)]);
% end
% min_max=(distance==min(distance)) | (distance==max(distance));
measure_mirror_real=[x_distance y_distance];
dimension_pattern=[floor(x_distance/5) floor(y_distance/5)]*5;
delta=floor([measure_mirror_real-dimension_pattern])/2;
%delta=[0 0];
n=floor([(dimension_pattern(1))/(pattern_size(2)-1) (dimension_pattern(2))/(pattern_size(1)-1)]);

n=mean(n);

5+delta(1):n:((pattern_size(2))*n);
5+delta(2):n:((pattern_size(1))*n);




xxp=[];
yyp=[];
cnt_xxyyp=0;
for xx=5+delta(1):n:((pattern_size(1))*n)
    for yy=5+delta(2):n:((pattern_size(2))*n)
        cnt_xxyyp=cnt_xxyyp+1;
        xxp(cnt_xxyyp)=xx;
        yyp(cnt_xxyyp)=yy;
    end
end
figure;
plot(xxp, yyp,'or')
hold on
y_distance_plus=floor(y_distance_plus/5)*5;
cnt_cut=0;
for xxs=(xx+n):n:y_distance_plus
    cnt_cut=cnt_cut+1;
    %aux=
    for yys=(5+delta(2)+n*cnt_cut):n:((pattern_size(2)*n)-n*cnt_cut)
        plot(xxs, yys,'or')
        cnt_xxyyp=cnt_xxyyp+1;
        xxp(cnt_xxyyp)=xxs;
        yyp(cnt_xxyyp)=yys;
%         xxs
%         yys
%        pause(0.5);
            %cnt_xxyyp=cnt_xxyyp+1;
    end
end

figure;
plot3(model_points{i}(1,:),model_points{i}(2,:),model_points{i}(3,:),".b")
hold on
plot3(xxp,yyp,zeros(size(yyp)),'or')
zlim([-1 1])
grid on
xlabel('x (mm)')
ylabel('y (mm)')
zlabel('z (mm)')
title(['Feature points on the object plane'])



%% Grid in model points for figures in paper
addpath('IR\')
figure;
points_c={};
points_c_no_dist={};
points_p={};
corners_model={};
for i=1:n_poses
    
    disp(['Pose ' num2str(i-1)])
    obj_x=xxp;
    obj_y=yyp;
%     ps=search_x_axis(model_points{i}(1,:)',model_points{i}(2,:)',model_points{i}(3,:)');
%     pd=search_y_axis(model_points{i}(1,:)',model_points{i}(2,:)',model_points{i}(3,:)');
%     corners_3d=[ps pd];

% Buscando esquinas en cada pose
    corners_c1_aux=pinv(K1)*[double(corner_store{i})'; ones(1,size(corner_store{i}',2))];
    corners_c2_aux=pinv(K2)*[double(corner_store2{i})'; ones(1,size(corner_store{i}',2))];
    corners_c1_aux=corners_c1_aux(1:2,:);
    corners_c2_aux=corners_c2_aux(1:2,:);
    corners_c1_undist =comp_distortion_oulu(corners_c1_aux,K1_dist);
    corners_c2_undist =comp_distortion_oulu(corners_c2_aux,K2_dist);
    corners_c1_undist=K1*[corners_c1_undist; ones(1,size(corners_c1_undist,2))];
    corners_c2_undist=K2*[corners_c2_undist; ones(1,size(corners_c2_undist,2))];
    corners_c1_undist=corners_c1_undist(1:2,:);
    corners_c2_undist=corners_c2_undist(1:2,:);
    data=pyrunfile("triangulate.py","X",K1=py.numpy.array(K1), K2=py.numpy.array(K2), R=py.numpy.array(R_st), t=py.numpy.array(t_st),pk1=py.numpy.array(corners_c1_undist),pk2=py.numpy.array(corners_c2_undist));
    corners_3d= double(data)';

%     corners_c1 = undistortPoints(double(corner_store{i}),cam1_params);
%     corners_c2 = undistortPoints(double(corner_store2{i}),cam2_params);
%     corners_3d = triangulate(corners_c1,corners_c2,StereoParams);
    B=planes(i,:);
    
    C = B(1)*corners_3d(:,1) + B(2)*corners_3d(:,2)  + B(3)*ones(size(corners_3d(:,1)));
    corners_3d_ideal=[corners_3d(:,1) corners_3d(:,2) C];
%     figure(2);
%     plot3(corners_3d_ideal(:,1),corners_3d_ideal(:,2),corners_3d_ideal(:,3),'ok')
%     xlabel('x')
%     ylabel('y')
%     zlabel('z')
    
    corners_3d_model=H_stores{i}*(corners_3d_ideal'-corners_3d_ideal(1,:)');
    corners_model{i}=corners_3d_model;
%     figure;
%     plot3(corners_3d_model(1,:),corners_3d_model(2,:),corners_3d_model(3,:),'ok')
%     xlabel('x')
%     ylabel('y')
%     zlabel('z')
%     zlim([-1 1])

% Final de busqueda de esquinas
    subplot(221)
%    figure(1);
%corners_3d_ideal
    %corners_3d_model=corners_3d_model';
    corners_3d1=[corners_3d_model corners_3d_model(:,1)];
    plot3(corners_3d1(1,:),corners_3d1(2,:),corners_3d1(3,:))
    hold on
    %plot3(obj_x(1),obj_y(1),0,'*k')
    plot3(obj_x(1:end),obj_y(1:end),zeros(size(obj_y(1:end))),'or','MarkerSize',5)
    % plot3(obj_x(2:end),obj_y(2:end),zeros(12*16-1,1),'*r')
    f_points_object=[xxp; yyp; zeros(size(yyp))];
    % f_points_object=[x_pattern(:)'; y_pattern(:)'; zeros(12*16,1)'];
    hold off
    zlim([-1 1])
    grid on
    xlabel('x (mm)')
    ylabel('y (mm)')
    zlabel('z (mm)')
    view(0,90)
    %title('Feature points on the object plane')
    subplot(222)
    %figure;
    H_aux=pinv(H_stores{i});
    f_points_cam=H_aux*f_points_object+p0_store(:,i);
    corners_3d1_cam=H_aux*corners_3d1+p0_store(:,i);
    plot3(f_points_cam(1,:),f_points_cam(2,:),f_points_cam(3,:),'or','MarkerSize',5)
    hold on
    plot3(corners_3d1_cam(1,:),corners_3d1_cam(2,:),corners_3d1_cam(3,:),'b')
    hold off
    grid on
    %zlim([-1 1])
    xlabel('x (mm)')
    ylabel('y (mm)')
    zlabel('z (mm)')
    axis equal
    view([45 45])
    %title('Feature points in the camera coordinate system')
    sgtitle(['Pose ' num2str(i-1)])
%     b0=[0,0,0]';
%     bx=[11,0,0]';
%     by=[0,11,0]';
%     H_aux=pinv(H_stores{i});
%     f_points=H_aux*[b0 bx by]+p0_store(:,i);
%     trasl=f_points(:,1);
%     dist_x=f_points(:,2)-f_points(:,1);
%     dist_y=f_points(:,3)-f_points(:,1);
%     nx=dist_x/norm(dist_x);
%     ny=dist_y/norm(dist_y);
%     nz=cross(nx,ny);
%     R_m_c=[nx';ny';nz'];
%     %H_m_c=[R_m_c t'];
%     H_m_c=[R_m_c t];
    %points_c=K1*H_m_c*[f_points_cam; ones(1,12*16)]
    points_c_aux=K1*f_points_cam;
    points_c{i}=points_c_aux./points_c_aux(3,:);
%     if length(cam1_params.RadialDistortion)==2
%         k_dist=[cam1_params.RadialDistortion 0 cam1_params.TangentialDistortion];
%     else
%         k_dist=[cam1_params.RadialDistortion cam1_params.TangentialDistortion];     
%     end
    points_c_no_dist{i}=points_c{i};
    norm_pixel=pinv(K1)*points_c{i};
    coord=norm_pixel(1:2,:);
    points_c{i}=distortion_points(K1,K1_dist,coord);
    %sum(points_c{i}-points_c_no_dist{i})
    subplot(223);
    %figure;
%     if i==1
%         imshow(zeros(height,width))  
%     end
    bg_screen=ones(height,width);
    bg_screen(1:3,:)=0;
    bg_screen(end-2:end,:)=0;
    bg_screen(:,1:3)=0;
    bg_screen(:,end-2:end)=0;
    imshow(bg_screen) 
    hold on
    plot(points_c{i}(1,:),points_c{i}(2,:),'ob','MarkerSize',6)
    %hold off
    %axis equal
    %title('Feature points in camera')



%     [row, col]=find(mask_final2{i}==1);
%     %figure; imshow(mask_final2{i})
%     x2=col;
%     y2=row;
%     phi2x=[];
%     phi2y=[];
%     for m=1:length(x2)
%         phi2x(m)=f2x(y2(m),x2(m));
%         phi2y(m)=f2y(y2(m),x2(m));
%     end


        % Load point coordinates
    f1x=fase_c1{i,1};
    f1y=fase_c1{i,2};
    [row, col]=find(mask_final{i}==1);
    x=col;
    y=row;
    phi1x=[];
    phi1y=[];
    for m=1:length(x)
        phi1x(m)=f1x(y(m),x(m));
        phi1y(m)=f1y(y(m),x(m));
    end

    coord_x=points_c{i}(1,:);
    coord_y=points_c{i}(2,:);
    %Interpolamos fase en x
    f1x_int= griddata(double(x),double(y),phi1x,double(coord_x),double(coord_y),"cubic");
    f1y_int= griddata(double(x),double(y),phi1y,double(coord_x),double(coord_y),"cubic");
% ojo
%     [x,y] = meshgrid(1:width, 1:height);   
%     PhasexPoints = interp2(double(x),double(y),f1x,double(coord_x),double(coord_y));
%     PhaseyPoints = interp2(double(x),double(y),f1y,double(coord_x),double(coord_y));
%     format long
%     points_p{i}(:,1)=P*PhasexPoints/(2*pi);
%     points_p{i}(:,2)=P*PhaseyPoints/(2*pi);
%ojo
    format long
    points_p{i}(:,1)=P*f1x_int/(2*pi);
    points_p{i}(:,2)=P*f1y_int/(2*pi);
    subplot(224);
    %figure;
%     if i==1
%         imshow(zeros(H,W))
%     end
    bg_screen=ones(H,W);
    bg_screen(1:2,:)=0;
    bg_screen(end-1:end,:)=0;
    bg_screen(:,1:2)=0;
    bg_screen(:,end-1:end)=0;

    imshow(bg_screen)
    hold on
    plot(points_p{i}(:,1),points_p{i}(:,2),'ok','MarkerSize',6)
    %title('Feature points in projector')
    pause(2);
end



%%
i=32;
figure;
    disp(['Pose ' num2str(i-1)])
    obj_x=xxp;
    obj_y=yyp;
%     ps=search_x_axis(model_points{i}(1,:)',model_points{i}(2,:)',model_points{i}(3,:)');
%     pd=search_y_axis(model_points{i}(1,:)',model_points{i}(2,:)',model_points{i}(3,:)');
%     corners_3d=[ps pd];

% Buscando esquinas en cada pose
    corners_c1_aux=pinv(K1)*[double(corner_store{i})'; ones(1,size(corner_store{i}',2))];
    corners_c2_aux=pinv(K2)*[double(corner_store2{i})'; ones(1,size(corner_store{i}',2))];
    corners_c1_aux=corners_c1_aux(1:2,:);
    corners_c2_aux=corners_c2_aux(1:2,:);
    corners_c1_undist =comp_distortion_oulu(corners_c1_aux,K1_dist);
    corners_c2_undist =comp_distortion_oulu(corners_c2_aux,K2_dist);
    corners_c1_undist=K1*[corners_c1_undist; ones(1,size(corners_c1_undist,2))];
    corners_c2_undist=K2*[corners_c2_undist; ones(1,size(corners_c2_undist,2))];
    corners_c1_undist=corners_c1_undist(1:2,:);
    corners_c2_undist=corners_c2_undist(1:2,:);
    data=pyrunfile("triangulate.py","X",K1=py.numpy.array(K1), K2=py.numpy.array(K2), R=py.numpy.array(R_st), t=py.numpy.array(t_st),pk1=py.numpy.array(corners_c1_undist),pk2=py.numpy.array(corners_c2_undist));
    corners_3d= double(data)';

%     corners_c1 = undistortPoints(double(corner_store{i}),cam1_params);
%     corners_c2 = undistortPoints(double(corner_store2{i}),cam2_params);
%     corners_3d = triangulate(corners_c1,corners_c2,StereoParams);
    B=planes(i,:);
    
    C = B(1)*corners_3d(:,1) + B(2)*corners_3d(:,2)  + B(3)*ones(size(corners_3d(:,1)));
    corners_3d_ideal=[corners_3d(:,1) corners_3d(:,2) C];
%     figure(2);
%     plot3(corners_3d_ideal(:,1),corners_3d_ideal(:,2),corners_3d_ideal(:,3),'ok')
%     xlabel('x')
%     ylabel('y')
%     zlabel('z')
    
    corners_3d_model=H_stores{i}*(corners_3d_ideal'-corners_3d_ideal(1,:)');
    corners_model{i}=corners_3d_model;
%     figure;
%     plot3(corners_3d_model(1,:),corners_3d_model(2,:),corners_3d_model(3,:),'ok')
%     xlabel('x')
%     ylabel('y')
%     zlabel('z')
%     zlim([-1 1])

% Final de busqueda de esquinas
    subplot(221)
%    figure(1);
%corners_3d_ideal
    %corners_3d_model=corners_3d_model';
    corners_3d1=[corners_3d_model corners_3d_model(:,1)];
    plot3(corners_3d1(1,:),corners_3d1(2,:),corners_3d1(3,:))
    hold on
    %plot3(obj_x(1),obj_y(1),0,'*k')
    plot3(obj_x(1:end),obj_y(1:end),zeros(size(obj_y(1:end))),'or','MarkerSize',4)
    % plot3(obj_x(2:end),obj_y(2:end),zeros(12*16-1,1),'*r')
    f_points_object=[xxp; yyp; zeros(size(yyp))];
    % f_points_object=[x_pattern(:)'; y_pattern(:)'; zeros(12*16,1)'];
    hold off
    zlim([-1 1])
    grid on
    xlabel('x (mm)')
    ylabel('y (mm)')
    zlabel('z (mm)')
    view(45,45)
    %title('Feature points on the object plane')
    subplot(222)
    %figure;
    H_aux=pinv(H_stores{i});
    f_points_cam=H_aux*f_points_object+p0_store(:,i);
    corners_3d1_cam=H_aux*corners_3d1+p0_store(:,i);
    plot3(f_points_cam(1,:),f_points_cam(2,:),f_points_cam(3,:),'or','MarkerSize',4)
    hold on
    plot3(corners_3d1_cam(1,:),corners_3d1_cam(2,:),corners_3d1_cam(3,:),'b')
    hold off
    grid on
    %zlim([-1 1])
    xlabel('x (mm)')
    ylabel('y (mm)')
    zlabel('z (mm)')
    axis equal
    view([45 45])
    %title('Feature points in the camera coordinate system')
    sgtitle(['Pose ' num2str(i-1)])
%     b0=[0,0,0]';
%     bx=[11,0,0]';
%     by=[0,11,0]';
%     H_aux=pinv(H_stores{i});
%     f_points=H_aux*[b0 bx by]+p0_store(:,i);
%     trasl=f_points(:,1);
%     dist_x=f_points(:,2)-f_points(:,1);
%     dist_y=f_points(:,3)-f_points(:,1);
%     nx=dist_x/norm(dist_x);
%     ny=dist_y/norm(dist_y);
%     nz=cross(nx,ny);
%     R_m_c=[nx';ny';nz'];
%     %H_m_c=[R_m_c t'];
%     H_m_c=[R_m_c t];
    %points_c=K1*H_m_c*[f_points_cam; ones(1,12*16)]
    points_c_aux=K1*f_points_cam;
    points_c{i}=points_c_aux./points_c_aux(3,:);
%     if length(cam1_params.RadialDistortion)==2
%         k_dist=[cam1_params.RadialDistortion 0 cam1_params.TangentialDistortion];
%     else
%         k_dist=[cam1_params.RadialDistortion cam1_params.TangentialDistortion];     
%     end
    points_c_no_dist{i}=points_c{i};
    norm_pixel=pinv(K1)*points_c{i};
    coord=norm_pixel(1:2,:);
    points_c{i}=distortion_points(K1,K1_dist,coord);
    %sum(points_c{i}-points_c_no_dist{i})
    subplot(223);
    %figure;
%     if i==1
%         imshow(zeros(height,width))  
%     end
    bg_screen=ones(height,width);
    bg_screen(1:3,:)=0;
    bg_screen(end-2:end,:)=0;
    bg_screen(:,1:3)=0;
    bg_screen(:,end-2:end)=0;
    imshow(bg_screen) 
    hold on
    plot(points_c{i}(1,:),points_c{i}(2,:),'ob','MarkerSize',4)
    %hold off
    %axis equal
    %title('Feature points in camera')



%     [row, col]=find(mask_final2{i}==1);
%     %figure; imshow(mask_final2{i})
%     x2=col;
%     y2=row;
%     phi2x=[];
%     phi2y=[];
%     for m=1:length(x2)
%         phi2x(m)=f2x(y2(m),x2(m));
%         phi2y(m)=f2y(y2(m),x2(m));
%     end


        % Load point coordinates
    f1x=fase_c1{i,1};
    f1y=fase_c1{i,2};
    [row, col]=find(mask_final{i}==1);
    x=col;
    y=row;
    phi1x=[];
    phi1y=[];
    for m=1:length(x)
        phi1x(m)=f1x(y(m),x(m));
        phi1y(m)=f1y(y(m),x(m));
    end

    coord_x=points_c{i}(1,:);
    coord_y=points_c{i}(2,:);
    %Interpolamos fase en x
    f1x_int= griddata(double(x),double(y),phi1x,double(coord_x),double(coord_y),"cubic");
    f1y_int= griddata(double(x),double(y),phi1y,double(coord_x),double(coord_y),"cubic");
% ojo
%     [x,y] = meshgrid(1:width, 1:height);   
%     PhasexPoints = interp2(double(x),double(y),f1x,double(coord_x),double(coord_y));
%     PhaseyPoints = interp2(double(x),double(y),f1y,double(coord_x),double(coord_y));
%     format long
%     points_p{i}(:,1)=P*PhasexPoints/(2*pi);
%     points_p{i}(:,2)=P*PhaseyPoints/(2*pi);
%ojo
    format long
    points_p{i}(:,1)=P*f1x_int/(2*pi);
    points_p{i}(:,2)=P*f1y_int/(2*pi);
    subplot(224);
    %figure;
%     if i==1
%         imshow(zeros(H,W))
%     end
    bg_screen=ones(H,W);
    bg_screen(1:2,:)=0;
    bg_screen(end-1:end,:)=0;
    bg_screen(:,1:2)=0;
    bg_screen(:,end-1:end)=0;

    imshow(bg_screen)
    hold on
    plot(points_p{i}(:,1),points_p{i}(:,2),'ok','MarkerSize',4)
    %title('Feature points in projector')
    pause(0.5);
%% Cartagena
% Caracteristica patron
%patternSize = pattern_size;  %[12, 18];       
worldPoints = [obj_x' obj_y'];


M = zeros([size(worldPoints,1),2,n_poses,2]);

supr_poses=[];
correct_poses=[];
cnt_poses_erroneas=0;
cnt_poses_correctas=0;
for i=1:n_poses
    if sum(sum(isnan(points_p{i})))~=0
        cnt_poses_erroneas=cnt_poses_erroneas+1;
        supr_poses(cnt_poses_erroneas)=i;
    else
        cnt_poses_correctas=cnt_poses_correctas+1;
        correct_poses(cnt_poses_correctas)=i;
    end   
end
disp(['Non adequate poses: ' num2str(supr_poses)])
disp(['Adequate poses: ' num2str(correct_poses)])

for i = 1:n_poses
    disp(['Processing pose ' num2str(i)])
    aux_points_c=points_c{i};
    aux_points_p=points_p{i};
    M(:,:,i,1) = [aux_points_c(1,:)' aux_points_c(2,:)'];
    M(:,1,i,2) = aux_points_p(:,1);
    M(:,2,i,2) = aux_points_p(:,2);
end

%Estereo Calibration using Matlab Toolbox
NRcoef = 2;
SkewFlag = 0;
Flagtang=logical(1);
%[1 6 7 38 58 61 62 72 77 78 80]
correct_poses=[1 3:6 9:38 40:58 60:61 64:72 74:77 80];
correct_poses=[1 3:6 9:38 40:43 45:47 49:58 60:61 64:65 67:72 74:77 80];
[StereoParams1,imagesUsed,estimationErrors] = estimateCameraParameters(M(:,:,correct_poses,:),worldPoints,'EstimateSkew',logical(SkewFlag),'NumRadialDistortionCoefficients',NRcoef,'EstimateTangentialDistortion',Flagtang);
%[StereoParams,imagesUsed,estimationErrors] = estimateCameraParameters(M(:,:,[8:18 20:30],:),worldPoints,'EstimateSkew',logical(SkewFlag),'NumRadialDistortionCoefficients',NRcoef);

figure, showReprojectionErrors(StereoParams1);


res_K1=[width height];
res_K2=[H W];
K1
K1_dist
data=pyrunfile("stereocalibration_1cam_1fixed.py","data",M=py.numpy.array([M]), worldPoints=py.numpy.array([worldPoints]), correct_poses=py.numpy.array([correct_poses]), res_K1=py.numpy.array([res_K1]),res_K2=py.numpy.array([res_K2]),K1=py.numpy.array([K1]),K1_dist=py.numpy.array([K1_dist]));
    %res = pyrunfile("addac.py","z",x=3,y=2)
K1 = double(data{'K1'})
K1_dist = double(data{'K1_dist'})
Kproj= double(data{'K2'})
Kproj_dist = double(data{'K2_dist'})
R_proj = double(data{'R'})
t_proj= double(data{'t'})
E = double(data{'E'})
F= double(data{'F'})
rms2= double(data{'rms2'}) % Mean reprojection error K2
rett = double(data{'rett'}) % Mean reprojection error Extrinsic




%save([direc 'calib_sl_python_dig.mat'],"K1","K1_dist","Kproj","Kproj_dist","R_proj","t_proj","Kt","Kt_dist","R_th","t_th")
%save([direc 'calib_sl_python_dig.mat'],"K1","K1_dist","Kproj","Kproj_dist","R_proj","t_proj")





%%


%% VIS-IR Calibration
figure;
%corner_store_vis={};
for i=1:n_poses
%     %i=1;
%     if(i<11), aux = ['0' num2str(i-1)]; else, aux = num2str(i-1); end
     imshow(texture_img{i})
%     % hold on
%     % for j=1:size(corner_store{i},2)
%     %     plot(corner_store{i}(1,j),corner_store{i}(2,j),'.r')
%     % 
%     % end
%     Y=imsharpen(texture_img{i},'Radius',3,'Amount',2);
%     level = graythresh(Y);        %Umbralizamos
%     mask1 = imbinarize(Y,level);   
%     [mask_w_w,numWhite] = bwlabel(mask1);
%     mask_w_w = bwareaopen(mask_w_w,100000);  %Retirando manchas blancas
%     [L,numBlack] = bwlabel(~mask_w_w); %Etiquetamos los agrupamientos de pixeles.
%     mask= ~bwareaopen(L,100000);   %Retirando manchas negras
%     mirror_corners = pgonCorners(mask,4);
%    corner_store_vis{i}=mirror_corners;
 
    hold on
    plot(corner_store{i}(:,1),corner_store{i}(:,2),'+r')
    sgtitle(['Pose ' aux])
    pause(0.5);
end
close all
%%
addpath('C:\Users\Eberto Benjumea\Desktop\OneDrive - Universidad Tecnológica de Bolívar\Doctorado\Proyectos\Codigos Reconstruccion 3D Learning\Codigos Tesis All Digital Features\IR')

close all
img_ir=import_ir_images([direc IR],n_poses,";");
%img_ir=import_ir_images([direc IR],n_poses,"\t"); % import thermal images from .txt files
size_ir=size(img_ir{1});

corner_store_ir={};


%% Activo
close all
figure;
for i=1:n_poses
    if(i<11), aux = ['0' num2str(i-1)]; else, aux = num2str(i-1); end
    disp(['Pose ' num2str(i-1)])
%------------Processing IR images-----------------------------
    [mask_ir{i},corner_store_ir{i}]=find_mirror_nmask_ir(img_ir{i}, n_corners,1);
    pause(0.5);
end
close all

%% IR Initial calibration for distortion estimation
figure;
hold on;
corner1=[];
corner2=[];
corner3=[];
corner4=[];
corner5=[];
corner6=[];
for i=1:n_poses
    plot(corners_model{i}(1,:),corners_model{i}(2,:),'*')
    corner1(i,:)=corners_model{i}(1:2,1)';
    corner2(i,:)=corners_model{i}(1:2,2)';
    corner3(i,:)=corners_model{i}(1:2,3)';
    corner4(i,:)=corners_model{i}(1:2,4)';
    corner5(i,:)=corners_model{i}(1:2,5)';
    corner6(i,:)=corners_model{i}(1:2,6)';
end
median_corner1 = median(corner1); 
median_corner2 = median(corner2); 
median_corner3 = median(corner3); 
median_corner4 = median(corner4); 
median_corner5 = median(corner5); 
median_corner6 = median(corner6); 

worldpoints_thermal=[median_corner1; median_corner2; median_corner3;
                     median_corner4; median_corner5; median_corner6];


M=[];
for i=1:n_poses
    M(:,:,i)=corner_store_ir{i};
end
%IR Calibration using Matlab Toolbox

NRcoef = 2;
SkewFlag = 0;
Flagtang=logical(1);
correct_poses=[2:5 8:33 39:45 48:55 58 63:67 69:71 74:81];


%[params,imagesUsed,estimationErrors] = estimateCameraParameters(M(:,:,[1:16 18:27 29:34],:),worldPoints,'EstimateSkew',logical(SkewFlag),'NumRadialDistortionCoefficients',NRcoef,'EstimateTangentialDistortion',Flagtang);
[params_ir,imagesUsed,estimationErrors] = estimateCameraParameters(M(:,:,correct_poses),worldpoints_thermal,'EstimateSkew',logical(SkewFlag),'NumRadialDistortionCoefficients',NRcoef,'EstimateTangentialDistortion',Flagtang);
% poses_used=[2:5 8:33 38:45 47:55 57:58 60 63:67 69:71 74:81];
% poses_used(imagesUsed)
% [1:34 39:46 49:56 58:59 61 63:67 69:72 74:81]

figure, showReprojectionErrors(params_ir);
% A=params_ir;
% A.IntrinsicMatrix'
% distir_1=A.RadialDistortion
% distir_2=A.TangentialDistortion
% [distir_1 distir_2 0]

res_K2=size_ir ;
data=pyrunfile("calibration_1cam.py","data",M=py.numpy.array([M]), worldPoints=py.numpy.array([worldpoints_thermal]), correct_poses=py.numpy.array([correct_poses]),res_K2=py.numpy.array([res_K2]));
    %res = pyrunfile("addac.py","z",x=3,y=2)

Kt= double(data{'K2'})
Kt_dist = double(data{'K2_dist'})
rms2= double(data{'rms2'}) % Mean reprojection error Kt



%%

    i=32
    figure;
    imagesc(img_ir{i})
    hold on
    plot(corner_store_ir{i}(:,1),corner_store_ir{i}(:,2),'+k','MarkerSize',12)




%% Test corners

close all
%corner_store_vis=corner_store;
homography_m={};
ir_points={};
figure;
for i=1:n_poses
    sgtitle(['Pose ' num2str(i-1)])
    subplot(121)
    imshow(texture_img{i})
    hold on
    subplot(122)
    imshow(img_ir{i},[])
    hold on

    for j=1:size(corner_store{i},1)
        subplot(121)
        plot(corner_store{i}(j,1),corner_store{i}(j,2),'+r')
        subplot(122)
        plot(corner_store_ir{i}(j,1),corner_store_ir{i}(j,2),'+r')
        
        %pause();
    end


    corners_c1_aux=pinv(K1)*[double(corner_store{i})'; ones(1,size(corner_store{i}',2))];
    corners_ir_aux=pinv(Kt)*[double(corner_store_ir{i})'; ones(1,size(corner_store_ir{i}',2))];
    corners_c1_aux=corners_c1_aux(1:2,:);
    corners_ir_aux=corners_ir_aux(1:2,:);
    corners_c1_undist =comp_distortion_oulu(corners_c1_aux,K1_dist);
    corners_ir_undist =comp_distortion_oulu(corners_ir_aux,Kt_dist);
    corners_c1_undist=K1*[corners_c1_undist; ones(1,size(corners_c1_undist,2))];
    corners_ir_undist=Kt*[corners_ir_undist; ones(1,size(corners_ir_undist,2))];
    corners_c1=corners_c1_undist(1:2,:)';
    corners_ir=corners_ir_undist(1:2,:)';


    %corners_c1 = undistortPoints(double(corner_store{i}),cam1_params);
    %corners_ir = undistortPoints(double(corner_store_ir{i}),params_ir);
    
    pin=[corners_c1(:,1)'; corners_c1(:,2)'];
    pout=[corners_ir(:,1)'; corners_ir(:,2)'];
    homography_m{i}= homography_solve(pin, pout);
    ir_points{i} = homography_transform(points_c_no_dist{i}(1:2,:),homography_m{i});
    %sum((vecnorm(pout-ir_points{i})))/6
%     if length(params_ir.RadialDistortion)==2
%         k_dist=[params_ir.RadialDistortion params_ir.TangentialDistortion 0];
%     else
%         k_dist=[params_ir.RadialDistortion params_ir.TangentialDistortion];     
%     end
    aux=ir_points{i};
    norm_pixel=pinv(Kt)*[aux; ones(1,size(aux,2))];
    coord=norm_pixel(1:2,:);
    ir_points{i} =distortion_points(Kt,Kt_dist,coord);
    pause(0.1);
end
close all
figure;
for i=1:n_poses
    subplot(121)
    imshow(texture_img{i})
    hold on
    plot(points_c{i}(1,:),points_c{i}(2,:),'+r')
    subplot(122)
    imshow(img_ir{i},[])
    hold on
    plot(ir_points{i}(1,:),ir_points{i}(2,:),'+r')
    sgtitle(['Pose ' num2str(i-1)])
    pause(0.5);
end
close all

%% Corners Homography without distortions
% close all
% %corner_store_vis=corner_store;
% homography_m={};
% ir_points={};
% figure;
% for i=1:n_poses
%     sgtitle(['Pose ' num2str(i-1)])
%     subplot(121)
%     imshow(texture_img{i})
%     hold on
%     subplot(122)
%     imshow(img_ir{i},[])
%     hold on
% 
%     for j=1:size(corner_store{i},1)
%         subplot(121)
%         plot(corner_store{i}(j,1),corner_store{i}(j,2),'+r')
%         subplot(122)
%         plot(corner_store_ir{i}(j,1),corner_store_ir{i}(j,2),'+r')
%         
%         %pause();
%     end
% 
%     pin=[corner_store{i}(:,1)'; corner_store{i}(:,2)'];
%     pout=[corner_store_ir{i}(:,1)'; corner_store_ir{i}(:,2)'];
%     homography_m{i}= homography_solve(pin, pout);
%     ir_points{i} = homography_transform(points_c{i}(1:2,:),homography_m{i});
%     pause();
% end
% close all
% figure;
% for i=1:n_poses
%     subplot(121)
%     imshow(texture_img{i})
%     hold on
%     plot(points_c{i}(1,:),points_c{i}(2,:),'+r')
%     subplot(122)
%     imshow(img_ir{i},[])
%     hold on
%     plot(ir_points{i}(1,:),ir_points{i}(2,:),'+r')
%     sgtitle(['Pose ' num2str(i-1)])
%     pause();
% end
% close all

%% Figure for paper
figure;
i=32;
subplot(121)
imshow(texture_img{i})
hold on
plot(points_c{i}(1,:),points_c{i}(2,:),'or')
plot(corner_store{i}(:,1), corner_store{i}(:,2),'+b','MarkerSize',12)
subplot(122)
imagesc(img_ir{i})
axis equal
ylim([1 192])
c = colorbar;
c.Label.String = '° C';

% imshow(img_ir{i},[])
% c = colorbar;
% c.Label.String = '° C';

hold on
plot(ir_points{i}(1,:),ir_points{i}(2,:),'or')
plot(corner_store_ir{i}(:,1), corner_store_ir{i}(:,2),'+k','MarkerSize',12)
%sgtitle(['Pose ' num2str(i-1)])
%% figure for paper 2
bg_screen=ones(192,256);
bg_screen(1:2,:)=0;
bg_screen(end-1:end,:)=0;
bg_screen(:,1:2)=0;
bg_screen(:,end-1:end)=0;
figure;
imshow(bg_screen,[])
hold on
plot(ir_points{i}(1,:),ir_points{i}(2,:),'or')
plot(corner_store_ir{i}(:,1), corner_store_ir{i}(:,2),'+b','MarkerSize',12)



%%
numbers=[3:5 8:33 39:45 49:55 58 64:67 69:71 74:81];
for j=1:length(numbers)
    i=numbers(j);
    %i=i-1;
    disp(['Pose ' num2str(i)])

    subplot(121)
    imshow(texture_img{i})
    hold on
    plot(points_c{i}(1,:),points_c{i}(2,:),'+r')
    subplot(122)
    imshow(img_ir{i},[])
    hold on
    plot(ir_points{i}(1,:),ir_points{i}(2,:),'+r')
    sgtitle(['Pose ' num2str(i-1)])
    pause(0.1);
end


%%
M=[];
for i=1:n_poses
    M(:,:,i,1)=points_c{i}(1:2,:)';
    M(:,:,i,2)=ir_points{i}(1:2,:)';
end
%VIS-IR Calibration using Matlab Toolbox
correct_poses=[3:5 8:33 39:45 49:55 58 64:67 69:71 74:81];
%save([direc 'data_visir_dist_calib.mat'],'K','K_dist','Kt','kt_dist','M','correct_poses','worldPoints')
NRcoef = 2;
SkewFlag = 0;
Flagtang=logical(1);

%[params,imagesUsed,estimationErrors] = estimateCameraParameters(M(:,:,[1:16 18:27 29:34],:),worldPoints,'EstimateSkew',logical(SkewFlag),'NumRadialDistortionCoefficients',NRcoef,'EstimateTangentialDistortion',Flagtang);
[params,imagesUsed,estimationErrors] = estimateCameraParameters(M(:,:,[3:5 9:33 39:45 49:55 58 64:67 69:71 74:81],:),worldPoints,'EstimateSkew',logical(SkewFlag),'NumRadialDistortionCoefficients',NRcoef,'EstimateTangentialDistortion',Flagtang);
% poses_used=[2:5 8:33 38:45 47:55 57:58 60 63:67 69:71 74:81];
% poses_used(imagesUsed)
% [1:34 39:46 49:56 58:59 61 63:67 69:72 74:81]

figure, showReprojectionErrors(params);


% visir=params;
% K1=visir.CameraParameters1.IntrinsicMatrix'
% dist1=visir.CameraParameters1.RadialDistortion
% dist12=visir.CameraParameters1.TangentialDistortion
% Kt=visir.CameraParameters2.IntrinsicMatrix'
% distt=visir.CameraParameters2.RadialDistortion
% distt2=visir.CameraParameters2.TangentialDistortion
% R=visir.RotationOfCamera2
% t=visir.TranslationOfCamera2

data=pyrunfile("stereocalibration.py","data",M=py.numpy.array([M]), worldPoints=py.numpy.array([worldPoints]), correct_poses=py.numpy.array([correct_poses]),K1=py.numpy.array([K1]), K1_dist=py.numpy.array([K1_dist]),K2=py.numpy.array([Kt]),K2_dist=py.numpy.array([Kt_dist]));
    %res = pyrunfile("addac.py","z",x=3,y=2)
length(correct_poses)
R_th= double(data{'R'})
t_th = double(data{'t'})
rett= double(data{'rett'}) % Mean reprojection error 




%save([direc 'calib_vis_ir_all_digital.mat'],"K1","dist1","Kt","distt","R","t")
%save([direc,'Comp_Calib_stereoparams_vis_ir_all_digitalv9',num2str(SkewFlag),'_',num2str(NRcoef),'_',num2str(Flagtang),'.mat'],'params')
correct_poses=[3:5 9:33 39:45 49:55 58 64:67 69:71 74:81];
%save([direc, 'data_visir_calib.mat'],'M','worldPoints','K','correct_poses','K_dist')

%%



K1
K1_dist
Kproj
Kproj_dist
R_proj
t_proj


Kt
Kt_dist
R_th
t_th

save([direc 'calib_all_python_dig2.mat'],'K1','K1_dist','Kproj','Kproj_dist','R_proj','t_proj','Kt','Kt_dist','R_th','t_th')
%save('march05')
%%
res_K1
res_K2=[size_ir(1) size_ir(2)];
data=pyrunfile("stereocalibration_1cam_1fixed.py","data",M=py.numpy.array([M]), worldPoints=py.numpy.array([worldPoints]), correct_poses=py.numpy.array([correct_poses]), res_K1=py.numpy.array([res_K1]),res_K2=py.numpy.array([res_K2]),K1=py.numpy.array([K1]),K1_dist=py.numpy.array([K1_dist]));
    %res = pyrunfile("addac.py","z",x=3,y=2)
K1 = double(data{'K1'})
K1_dist = double(data{'K1_dist'})
Kt= double(data{'K2'})
Kt_dist = double(data{'K2_dist'})
R_th = double(data{'R'})
t_th= double(data{'t'})
E = double(data{'E'})
F= double(data{'F'})
rms3= double(data{'rms2'}) % Mean reprojection error K2
rett = double(data{'rett'}) % Mean reprojection error Extrinsic

Kproj
Kproj_dist
R_proj
t_proj

%%

function v = homography_solve(pin, pout)
% HOMOGRAPHY_SOLVE finds a homography from point pairs
%   V = HOMOGRAPHY_SOLVE(PIN, POUT) takes a 2xN matrix of input vectors and
%   a 2xN matrix of output vectors, and returns the homogeneous
%   transformation matrix that maps the inputs to the outputs, to some
%   approximation if there is noise.
%
%   This uses the SVD method of
%   http://www.robots.ox.ac.uk/%7Evgg/presentations/bmvc97/criminispaper/node3.html
% David Young, University of Sussex, February 2008
if ~isequal(size(pin), size(pout))
    error('Points matrices different sizes');
end
if size(pin, 1) ~= 2
    error('Points matrices must have two rows');
end
n = size(pin, 2);
if n < 4
    error('Need at least 4 matching points');
end
% Solve equations using SVD
x = pout(1, :); y = pout(2,:); X = pin(1,:); Y = pin(2,:);
rows0 = zeros(3, n);
rowsXY = -[X; Y; ones(1,n)];
hx = [rowsXY; rows0; x.*X; x.*Y; x];
hy = [rows0; rowsXY; y.*X; y.*Y; y];
h = [hx hy];
if n == 4
    [U, ~, ~] = svd(h);
else
    [U, ~, ~] = svd(h, 'econ');
end
v = (reshape(U(:,9), 3, 3)).';
end


function y = homography_transform(x, v)
% HOMOGRAPHY_TRANSFORM applies homographic transform to vectors
%   Y = HOMOGRAPHY_TRANSFORM(X, V) takes a 2xN matrix, each column of which
%   gives the position of a point in a plane. It returns a 2xN matrix whose
%   columns are the input vectors transformed according to the homography
%   V, represented as a 3x3 homogeneous matrix.
q = v * [x; ones(1, size(x,2))];
%q = v * [x; y; ones(1, size(x,2))];
p = q(3,:);
y = [q(1,:)./p; q(2,:)./p];
end