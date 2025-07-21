clc
clearvars
close all
method=input('Which method do you prefer?\n Press 1 for Two Stages UTB\n Press 2 for Two Stages 2 CITEDI\n Press 3 for Digital Features CITEDI\n Press 4 for Digital Features 1 UTB \n Press 5 for Digital Features 2 UTB \n Your answer: ');

switch method
    case 1
        % Two Stages
        
        load('D:\R3D\Adquisiciones\Experimento 50\Two stage\SL2\SL2_TSC_citediMay021.mat') % Load structured light system parameters
        
        %load('D:\R3D\Adquisiciones\Experimento 50\Two stage\Calib_stereoparams_vis_ir_0_3_1.mat') % Load structured light system parameters
        load('D:\R3D\Adquisiciones\Experimento 50\Two stage\VIS-IR\MayCalib_stereoparams_vis_ir_1_3_1.mat') % Load VIS-IR system parameters
        disp('Two Stages')
    case 2
        % Two Stages 2
        %load('D:\R3D\Adquisiciones\Experimento 53\Two stages\SL\SL_TSC_citediMay031') % Load structured light system parameters
        %load('D:\R3D\Adquisiciones\Experimento 53\Two stages\VIS-IR\MayCalib_stereoparams_vis_ir_1_3_1.mat') % Load VIS-IR system parameters
        load('D:\R3D\Adquisiciones\Experimento 54\Two stages\calib_all_python_trad.mat')
            
        disp('Two Stages 2')
    case 3
        % Digital Features CITEDI
        load('D:\R3D\2023\multimodal\Experimento 44\Digital Features\Calib_TSC_0_3.mat')
        load('D:\R3D\2023\multimodal\Experimento 44\Digital Features\Calib_stereoparams_vis_ir_all_digital1_3_1.mat')
        disp('Digital Features CITEDI')
    case 4
        % Digital Features I UTB
        %load('D:\R3D\Adquisiciones\Experimento 53\Digital features\Comp_Calib_stereoparams_vis_ir_all_digitalv91_3_1.mat')
        %load('D:\R3D\Adquisiciones\Experimento 53\Digital features\Calib_TSCv9_0_3.mat')
        %load('D:\R3D\Adquisiciones\Experimento 54\Digital features\Comp_Calib_stereoparams_vis_ir_all_digitalv90_2_1.mat')
        load('D:\R3D\Adquisiciones\Experimento 54\Digital features\calib_all_python_dig.mat')
        
        
        disp('Digital Features I UTB')
    case 5
        % Digital Features II UTB
        load('D:\R3D\Adquisiciones\Experimento 54\Digital features\calib_all_python_dig.mat')
        disp('Digital Features II UTB')
    otherwise
        disp('Run again the code')
end




%D:\R3D\Adquisiciones\Experimento 45\Objeto\Objeto final
direc='D:\R3D\Adquisiciones\Experimento 45\Objeto\Objeto final\'; % Load objects
direc='D:\R3D\Adquisiciones\Experimento 45\Objeto 2\';
direc='D:\R3D\Adquisiciones\Experimento 53\Objects\';
%direc='D:\R3D\Adquisiciones\Experimento 53\Objects\';
direc='D:\R3D\Adquisiciones\Experimento 54\Objects\'
%direc='D:\R3D\Adquisiciones\Experimento 50\Objects\'
%direc='D:\R3D\Adquisiciones\Experimento 54\Objects\';
%direc='D:\R3D\Adquisiciones\Experimento 54\Map_cal';
VIS='SL\';
IR='IR\';
inicio=9;
n_poses=29;
name_texture='P18_F26';
format='.PNG';


%% Graficamos camaras vistas desde el sistema coordenado de la camara visible
% Color axis: X: red, Y: green, Z: blue.
% pose de c1
close all
M_c1=[1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
% pose de t
M_t=pinv([R_th t_th; 0 0 0 1]);
% pose de proyector
Rp=R_proj;
tp=t_proj;
M_p=pinv([Rp tp; 0 0 0 1]);
figure;
hold on
graph_cam(M_c1,1)
graph_cam(M_t,2)
graph_cam(M_p,3)
hold off
axis equal
grid on
xlabel('x (mm)')
ylabel('y (mm)')
zlabel('z (mm)')
title('Locations of the cameras')


disp(['Locations of the cameras respect to the VIS-Cam (x y z in mm):'])
disp([''])
disp(['VIS-Cam: ' num2str((M_c1*[0;0;0;1])')])
disp(['IR-Cam: ' num2str((M_t*[0;0;0;1])')])
disp(['Projector: ' num2str((M_p*[0;0;0;1])')])

disp(['Distances of the cameras respect to the VIS-Cam (mm):'])
disp([''])
v_c1=M_c1*[0;0;0;1];
v_p=M_p*[0;0;0;1];
v_t=M_t*[0;0;0;1];
disp(['VIS-Cam: ' num2str(norm(v_c1(1:3,:)))])
disp(['IR-Cam: ' num2str(norm(v_t(1:3,:)))])
disp(['Projector: ' num2str(norm(v_p(1:3,:)))])
pause(1);

%%
% R_c1=[1 0 0; 0 1 0; 0 0 1];
% t_c1=[0, 0, 0];
% pose_c1 = rigid3d(R_c1,t_c1);
% 
% figure;
% cam_vis = plotCamera('AbsolutePose',pose_c1,'Opacity',0,'Label','VIS-Cam','Size',100)
% 
% pose_t = rigid3d(R,t');
% 
% hold on
% cam_ir = plotCamera('AbsolutePose',pose_t,'Opacity',0,'Label','IR-Cam','Size',100)
% 
% pose_p= rigid3d(Rp,tp');
% 
% hold on
% cam_ir = plotCamera('AbsolutePose',pose_p,'Opacity',0,'Label','Proj','Size',100)
% axis equal
%% Objects to reconstruct
%clc

disp('Objects to reconstruct.')
disp(' ')

figure;
texture_img={};
disp('Pose: ')
if isprime(n_poses)==1
    f = factor(n_poses+1);
    f=f(1);
    c=(n_poses+1)/f;
else
    f = factor(n_poses);
    f=f(1);
    c=n_poses/f;
end
for i=inicio:n_poses
    if i<11
        aux=['0' num2str(i-1)];
    else
        aux=num2str(i-1);
    end
    disp(aux)
    % Leemos la imagen
    texture_img{i}=imread([direc VIS 'pose_' aux '\H\' name_texture format]);    %Accedemos a imagen de textura
    img=texture_img{i};
    width=size(img,2);
    height=size(img,1);
    subplot(f,c,i)
    imshow(imadjust(img))
    title(['Object ' aux])
    % Buscamos las esquinas y la mascara via umbralziación
    sgtitle('Objects ')
    pause(0.5)
    
end
%close all

%% Obtención de fase
disp('Phase obtaining.')
disp(' ')
mode_phase=input(['If you want to obtain the phase from the images, please press 1\n' ...
    'If you want to read the phase from previously processed .mat files, \nPress 2\n Your answer: ']);
if mode_phase==1
    disp('Obtaining phase from images...')
    NStep = 18;
    NBits = 7;
    carp=VIS;
    [fase_c1]=estimate_phases(direc,carp,NStep,NBits,n_poses,width,height,format,inicio);
    %carp=camera2;
    %fase_c2=estimate_phases(direc,carp,NStep,NBits,n_poses,width,height,format);
    disp('Finished!')
elseif mode_phase==2
    disp('Reading phase from .mat files...')
    fase_c1={};

    %fase_c2={};
    for i=inicio:n_poses
        if i<11, aux=['0' num2str(i-1)]; else, aux=num2str(i-1); end
        load([direc VIS VIS(1:end-1) '_pose_' aux '_phases.mat'])
        fase_c1{i,1}=Fx;
        fase_c1{i,2}=Fy;
        %load([direc camera2 '\camera 2_pose_' aux '_phases.mat'])
        %fase_c2{i,1}=Fx;
        %fase_c2{i,2}=Fy;
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
for i=inicio:n_poses
%     fase_c1{i,1}(~mask_final{i})=NaN;
%     fase_c1{i,2}(~mask_final{i})=NaN;
%     fase_c2{i,1}(~mask_final2{i})=NaN;
%     fase_c2{i,2}(~mask_final2{i})=NaN;
    title(['Pose ' num2str(i-1)])
    subplot(121)
    imagesc(fase_c1{i,1})
    axis equal
    subplot(122)
    imagesc(fase_c1{i,2})
    axis equal
%     subplot(223)
%     imagesc(fase_c2{i,1})
%     axis equal
%     subplot(224)
%     imagesc(fase_c2{i,2})
%     axis equal
    sgtitle(['Pose ' num2str(i-1)])
    pause(1);
end


%% 3D data obtaining Eberto's Version
close all
img_ir=import_ir_images([direc IR],n_poses,";"); % import thermal images from .txt files

addpath('C:\Users\Eberto Benjumea\Desktop\OneDrive - Universidad Tecnológica de Bolívar\Doctorado\Proyectos\Codigos Reconstruccion 3D Learning\Codes Raul Purdue\StereoCalibCodes');

%addpath('C:\Users\Beatriz\OneDrive - Universidad Tecnológica de Bolívar\Doctorado\Proyectos\Codigos Reconstruccion 3D Learning\Codes Raul Purdue\StereoCalibCodes')
Direc=direc;  %Directorio de fases objeto
file=VIS;
carp=VIS;

Kt
Kt_dist
% Kt=params.CameraParameters2.IntrinsicMatrix';
% radial=params.CameraParameters2.RadialDistortion;
% if length(radial)==3
%             kt_dist=[radial(1) radial(2) params.CameraParameters2.TangentialDistortion radial(3) ];
% elseif length(radial)==2
%             kt_dist=[radial(1) radial(2) params.CameraParameters2.TangentialDistortion];
% end

n_r3d=n_poses;
%n_r3d=15;

%Mask=ones(height,width);
%Mask=maskC;


read_mask=0;
suave=0;
flag_rgb=1;
% st=StereoParams;
% cam1_params=st.CameraParameters1;
% cam2_params=st.CameraParameters2;
P=18;
for i=1:n_poses
    subplot(f,c,i)
    imshow(img_ir{i},[])
end
for i=inicio-1:n_r3d-1
    disp(['Pose ' num2str(i)])
    %carp='camera 1\'; %Ojo


    if i<10
        name_phase=[Direc carp carp(:,1:end-1) '_pose_0' num2str(i) '_phases.mat'];
        if read_mask==0
        %Mask=imread([Direc carp 'pose_0' num2str(i) '\im_51.png']);
            Mask=imread([Direc carp 'pose_0' num2str(i) '\H' '\P18_F26.PNG']); % descomentar para archivos de raul
        else

            
            % Verificar si el archivo existe
            if exist([Direc carp 'mask_pose_0' num2str(i) '.png'], 'file') == 2
                disp('El archivo existe.');
                Mask=imread([Direc carp 'mask_pose_0' num2str(i) '.png']);
            else
                disp('El archivo no existe.');
                Mask=zeros(height,width);
            end
            
        end
    else
        name_phase=[Direc carp carp(:,1:end-1) '_pose_' num2str(i) '_phases.mat'];
        if read_mask==0
        %Mask=imread([Direc carp 'pose_' num2str(i) '\im_51.png']);
            Mask=imread([Direc carp 'pose_' num2str(i) '\H' '\P18_F26.PNG']); % descomentar para archivos de raul
        
        else
                        % Verificar si el archivo existe
            if exist([Direc carp 'mask_pose_' num2str(i) '.png'], 'file') == 2
                disp('El archivo existe.');
                Mask=imread([Direc carp 'mask_pose_' num2str(i) '.png']);
            else
                disp('El archivo no existe.');
                Mask=zeros(height,width);
            end
            
        end
    end



    load(name_phase)



    if read_mask==0
        B = imbinarize(Mask, graythresh(Mask));
        SE = strel("disk",5);
        B = imdilate(B,SE);
        B = imclose(B,SE);
        B = imclose(B,SE);
        % SE = strel("disk",15);
        Mask = imclose(B,SE);
%         Mask = imfill(Mask,"holes"); % Ojo




%         X = imadjust(Mask);
%         BW = X > 21;
%         radius = 3;
%         decomposition = 0;
%         se = strel('disk', radius, decomposition);
%         BW = imopen(BW, se);
%         Mask = imfill(BW, 'holes');
    end
    fx=Fx;
    fy=Fy;
    FringePitch=18;
    typeCalib=4;
    Fd='Fy';
%     [cx,cy]=meshgrid(1:width,1:height);
%     cx=cx(:);
%     cy=cy(:);
    %[Mask,row,col] = segmentImage(img_ir{i+1});
    %figure;
    %imshow(Mask)
    K1
    K1_dist
%     radial=StereoParams.CameraParameters1.RadialDistortion;
%     if length(radial)==3
%         k_dist=[radial(1) radial(2) StereoParams.CameraParameters1.TangentialDistortion radial(3) ];
%     elseif length(radial)==2
%         k_dist=[radial(1) radial(2) StereoParams.CameraParameters1.TangentialDistortion 0];
%     end
    if i<11, aux=['0' num2str(i)]; else, aux=num2str(i); end
    % Mask=imread([direc VIS 'mask' aux '.png']);
    %Mask=ones(height,width);
    [row, col]=find(Mask==1);
    x=col;
    y=row;
    Fx_aux=[];
    Fy_aux=[];
    for m=1:length(row)
        Fx_aux(m)=fx(y(m),x(m));
        Fy_aux(m)=fy(y(m),x(m));
    end
    %disp('Hello')


    %fx=fx(:);
    %fy=fy(:);
    
    Fx_aux = fx(Mask(:)>0);
    Fy_aux = fy(Mask(:)>0);
    px = P*Fx_aux/(2*pi);
    py = P*Fy_aux/(2*pi);
    coor_p=[px py];
%     figure; imshow(zeros(800,1280))
%     hold on
%     plot(coor_p(:,1),coor_p(:,2),'.r')
%     figure; imshow(zeros(1200,1920))
%     hold on
%     plot(coor_c1(:,1),coor_c1(:,2),'.r')
    c1_points=double([col row]);
    p1_points=double([px py]);
    % Ojo corregir normalizacion de puntos


    c1_points_aux=pinv(K1)*[c1_points'; ones(1,size(c1_points',2))];
    p1_points_aux=pinv(Kproj)*[p1_points'; ones(1,size(p1_points',2))];
    c1_points_aux=c1_points_aux(1:2,:);
    p1_points_aux=p1_points_aux(1:2,:);



    c1_points_undist = comp_distortion_oulu(c1_points',K1_dist);
    p1_points_undist = comp_distortion_oulu(p1_points',K1_dist);
    
    %coor_c1_undist=undistortPoints(c1_points,cam1_params);
    %coor_c2_undist=undistortPoints(p1_points,cam2_params);
    tic
    %Points = triangulate(c1_points_undist',p1_points_undist',StereoParams);
    %Mask=imread('D:\R3D\Adquisiciones\Experimento 50\Objects\SL\mask_pose_24.png')
    %[XcM,YcM,ZcM,t1]=StereoReconstruction(fx,fy,Mask,FringePitch,StereoParams,typeCalib,Fd);
    [XcM,YcM,ZcM,t1]=StereoReconstruction_withcalibpython(fx,fy,Mask,FringePitch,K1,Kproj,K1_dist,Kproj_dist,R_proj,t_proj,typeCalib,Fd);
  
    
    %k=StereoParams.CameraParameters1.IntrinsicMatrix';
    K1
    %%% TEST
    oo=K1*[XcM(:)'; YcM(:)'; ZcM(:)'];
    oo=oo./oo(3,:);
    figure; imshow(texture_img{i+1})
    hold on
    plot(oo(1,1:50:end),oo(2,1:50:end),'.b')
    
    %%%%%
    ZcM=remove_peaks(ZcM);
        %ZcM = imgaussfilt(ZcM);
        
    ZcM=medfilt2(ZcM,[9 9]);
%     XcM1=Points(:,1);
%     YcM1=Points(:,2);
%     ZcM1=Points(:,3);
    XcM=XcM(:);
    YcM=YcM(:);
    ZcM=ZcM(:);    
%

%     p_mine=k*[1 0 0 0;0 1 0 0; 0 0 1 0]*[XcM1'; YcM1'; ZcM1'; ones(1,length(XcM1))];
%     p_raul=k*[1 0 0 0;0 1 0 0; 0 0 1 0]*[XcM'; YcM'; ZcM'; ones(1,length(XcM))];
%     p_mine=p_mine./p_mine(3,:);
%     p_raul=p_raul./p_raul(3,:);
%     figure;
%     plot3(XcM1,YcM1,ZcM1)
%     figure;
%     plot3(XcM,YcM,ZcM)
%     figure;
%     title('Mine')
%     imshow(zeros(2000,3000))
%     hold on
%     plot(p_mine(1,:),p_mine(2,:),'.')
%     figure;
%     title('Raul')
%     imshow(zeros(2000,3000))
%     hold on
%     plot(p_raul(1,:),p_raul(2,:),'.')
    
    
% load('D:\R3D\2023\multimodal\Experimento 44\Two Stages 2\VIS-IR\intrinsic_parameters_ir.mat')

    %i
    %obtain_texture_interpol4(Direc, file, Kt,kt_dist,StereoParams,params, XcM,YcM,ZcM,img_ir{i+1},i,height,width,method);
    
    obtain_texture_interpol6_withdatapython(Direc, file, Kt,Kt_dist,K1,K1_dist,R_th,t_th, XcM,YcM,ZcM,img_ir{i+1},texture_img{i+1},i,height,width,method);
    toc

%     Mask(1:400,:)=0;
%     Mask(:,1:400)=0;
%     Mask(end-400:end,:)=0;
%     Mask(:,end-400:end)=0;  
    %Mask=ones(1024,1280);
    [XcM,YcM,ZcM,t1]=StereoReconstruction_withcalibpython(fx,fy,Mask,FringePitch,K1,Kproj,K1_dist,Kproj_dist,R_proj,t_proj,typeCalib,Fd);
  
    
       ZcM=remove_peaks(ZcM);
        %ZcM = imgaussfilt(ZcM);
        
    ZcM=medfilt2(ZcM,[9 9]);
   %ZcM=-ZcM;
    %texture_img_ir=obtain_texture_interpol(params,XcM,YcM,ZcM,img_ir{i+1});
% 
    figure;
           s = surf(-XcM,YcM,-ZcM,'FaceColor', 'interp',...
                        'EdgeColor', 'none',...
                        'FaceLighting', 'phong');
    
        %view(-160,80)
            set(gca, 'DataAspectRatio', [1, 1, 1])
            %            view(45, 45);

            axis equal;

            camlight right
            axis on
            grid on
            xlabel('x (mm)')
            ylabel('y (mm)')
            zlabel('z (mm)')
            view(45,45)
            view(165,60)
            %set(gca, 'Zdir', 'reverse')
            %set(gca, 'Xdir', 'reverse')
            %zt = get(gca, 'ZTick');
            %set(gca, 'ZTick',zt, 'ZTickLabel',fliplr(zt))

            %set(gca, 'Zdir', 'reverse')
%    pause(5);
    close all
% 


end

%%
function graph_cam(pose_cam,type)
    axes =  [0, 0, 0, 1; 100, 0, 0, 1; 0, 100, 0, 1; 0, 0, 100, 1]';
    p2d = pose_cam(1:3, 1:4) * axes;
    line([p2d(1,1) p2d(1,2)],[p2d(2,1) p2d(2,2)],[p2d(3,1) p2d(3,2)],'Color','red');
    hold on
    line([p2d(1,1) p2d(1,3)],[p2d(2,1) p2d(2,3)],[p2d(3,1) p2d(3,3)],'Color','green');
    hold on
    line([p2d(1,1) p2d(1,4)],[p2d(2,1) p2d(2,4)],[p2d(3,1) p2d(3,4)],'Color','blue');
    if type==1
        text(p2d(1,1),p2d(2,1),p2d(3,1),'VIS-Cam')
    elseif type==2
        text(p2d(1,1),p2d(2,1),p2d(3,1),'IR-Cam')
    elseif type==3
        text(p2d(1,1),p2d(2,1),p2d(3,1),'Proj')
    end
end



