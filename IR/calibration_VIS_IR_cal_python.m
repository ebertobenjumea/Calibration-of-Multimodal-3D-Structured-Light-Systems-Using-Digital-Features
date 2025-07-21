close all
clear all
clc

addpath('C:\Users\Eberto Benjumea\Desktop\OneDrive - Universidad Tecnológica de Bolívar\Doctorado\Proyectos\Codigos Reconstruccion 3D Learning\Codigos Tesis All Digital Features\IR')
% direc='D:\R3D\2023\multimodal\Experimento 44\Two Stages\VIS-IR\'; % Two Stages 2
direc='D:\R3D\Adquisiciones\Experimento 54\Two stages\VIS-IR\';
vis='VIS\';
ir='IR\';
n_poses=85;


patternSize = [7, 9];%[5, 7];%[7, 9];  %[12, 18];         % [7, 21]
centerDistance = 30;%40; %30;            % 20
Pattern_type = 'asymmetric';  %'symmetric';     % asymmetric
Circ_color_vis = 'black';           % black
Circ_color_ir = 'white';           % white
img_ir=import_ir_images([direc ir],n_poses,";"); % import thermal images from .txt files
size_ir=size(img_ir{1});

M=[];
figure;
for i=1:n_poses
    if(i<11), aux = ['0' num2str(i-1)]; else, aux = num2str(i-1); end
    disp(['Pose ' num2str(i-1)])
%------------Processing VIS images-----------------------------
    img1=imread([direc vis 'pose_' aux '.PNG']);
    size_vis=size(img1);
    subplot(331)
    imshow(img1)
    title('VIS image')
    img=imadjust(img1);
    subplot(332)
    imshow(img)
    title('Adjusted VIS image')
    puntos_VIS = detectCircleGridPoints(img,patternSize,PatternType=Pattern_type,CircleColor=Circ_color_vis);
    disp(['Detected VIS points: ' num2str(size(puntos_VIS,1))])
    if size(puntos_VIS,1)==patternSize(1)*patternSize(2) 
        subplot(333) 
        imshow(img)
        hold on
        plot(puntos_VIS(:,1),puntos_VIS(:,2),'+r')
        title('Detected points on VIS image')
        M(:,:,i,1)=puntos_VIS;
    end
    
%------------Processing IR images-----------------------------

    subplot(334)
    imshow(img_ir{i},[])
    title('IR image')
    subplot(337)
    imhist(uint8(img_ir{i}))
    X=imadjust(uint8(img_ir{i}));
    subplot(335)
    imshow(X,[])
    title('Adjusted IR image')
    subplot(338)
    imhist(X)
    subplot(336)
    imshow(img_ir{i},[])
    hold on
    puntos_IR = detectCircleGridPoints(X,patternSize,PatternType=Pattern_type,CircleColor=Circ_color_ir);
    disp(['Detected IR points: ' num2str(size(puntos_IR,1))])
    if size(puntos_IR,1)==patternSize(1)*patternSize(2)  
        plot(puntos_IR(:,1),puntos_IR(:,2),'+r')
        title('Detected points on IR image')
        M(:,:,i,2)=puntos_IR;
    end
    hold off 

    sgtitle(['Pose ' aux])
    pause(0.5);
end

%% Calibration
worldPoints = generateCircleGridPoints(patternSize,centerDistance,'PatternType', Pattern_type);
NRcoef = 2;
SkewFlag = 0;
Flagtang=logical(1);
%save('data_calib_two_stages_calib2.mat',"M","worldPoints")

%% VIS calibration
A=M(:,:,:,1);
% params_vis= estimateCameraParameters(A(:,:,[1:4 6:13 17:18  23:25 28:32]),worldPoints,'EstimateSkew',logical(SkewFlag),'NumRadialDistortionCoefficients',NRcoef,'EstimateTangentialDistortion',Flagtang,'ImageSize',size_vis); %Two Stages 2
params_vis= estimateCameraParameters(A(:,:,:),worldPoints,'EstimateSkew',logical(SkewFlag),'NumRadialDistortionCoefficients',NRcoef,'EstimateTangentialDistortion',Flagtang,'ImageSize',size_vis);  % Two Stages
figure;
showReprojectionErrors(params_vis)
A=params_vis;
K1=A.IntrinsicMatrix'
distvis_1=A.RadialDistortion
distvis_2=A.TangentialDistortion

%% IR Calibration
B=M(:,:,:,2);
% params_ir= estimateCameraParameters(B(:,:,[1:4 6:7 9:13 16:21 23:32]),worldPoints,'EstimateSkew',logical(SkewFlag),'NumRadialDistortionCoefficients',NRcoef,'EstimateTangentialDistortion',Flagtang,'ImageSize',size_ir); % Two Stages 2
params_ir= estimateCameraParameters(B(:,:,:),worldPoints,'EstimateSkew',logical(SkewFlag),'NumRadialDistortionCoefficients',NRcoef,'EstimateTangentialDistortion',Flagtang,'ImageSize',size_ir); % Two Stages

figure;
showReprojectionErrors(params_ir)
B=params_ir;
Kt=B.IntrinsicMatrix'
distir_1=B.RadialDistortion
distir_2=B.TangentialDistortion

%save([direc 'intrinsic_parameters_ir.mat'],"Kt","distir_1","distir_2")

%% VIS-IR calibration
close all
SkewFlag = 0;
Flagtang=logical(1);
% params= estimateCameraParameters(M(:,:,[1:11
% 13:35],:),worldPoints,'EstimateSkew',logical(SkewFlag),'NumRadialDistortionCoefficients',NRcoef,'EstimateTangentialDistortion',Flagtang);
% Two Stages 2

[params,pairsUsed,estimationErrors]= estimateCameraParameters(M(:,:,:,:),worldPoints,'EstimateSkew',logical(SkewFlag),'NumRadialDistortionCoefficients',NRcoef,'EstimateTangentialDistortion',Flagtang);
clc
poses_calibrate=1:n_poses;
correct_poses=poses_calibrate(pairsUsed);
%poses_calibrate=poses_calibrate([1:17 20 24:length(poses_calibrate)]);
% 
[params,pairsUsed,estimationErrors]= estimateCameraParameters(M(:,:,poses_calibrate,:),worldPoints,'EstimateSkew',logical(SkewFlag),'NumRadialDistortionCoefficients',NRcoef,'EstimateTangentialDistortion',Flagtang);
% 
% list1=[1:17 19:20 22:25 27 29:30 32 34];
% list1(pairsUsed)
% [params,pairsUsed,estimationErrors]= estimateCameraParameters(M(:,:,[1:10 12:15 19 22:23 25 27 ],:),worldPoints,'EstimateSkew',logical(SkewFlag),'NumRadialDistortionCoefficients',NRcoef,'EstimateTangentialDistortion',Flagtang);
%%

res_K1=[size_vis(1) size_vis(2)];
res_K2=[size_ir(1) size_ir(2)];
load('D:\R3D\Adquisiciones\Experimento 54\Two stages\calib_sl_python_trad.mat')
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



%save([direc 'calib_all_python_trad.mat'],"K1","K1_dist","Kproj","Kproj_dist","R_proj","t_proj","Kt","Kt_dist","R_th","t_th")

%%
% SkewFlag = 0;
% [params,pairsUsed,estimationErrors]= estimateCameraParameters(M(:,:,[1:3 5:11 14:16 18:19 24 26:end],:),worldPoints,'EstimateSkew',logical(SkewFlag),'NumRadialDistortionCoefficients',NRcoef,'EstimateTangentialDistortion',Flagtang);
% list1=1:size(M(:,:,:,1),3);
% list1(pairsUsed)
% 
% figure;
% showReprojectionErrors(params)
% %save([direc,'JunCalib_stereoparams_vis_ir_',num2str(SkewFlag),'_',num2str(NRcoef),'_',num2str(Flagtang),'.mat'],'params')
% 


%%

% % pairs_used=[1:10 12:15 19 22:23 25 27 ];
% for j=1:length(poses_calibrate)
%     i=poses_calibrate(j);
%     if(i<11), aux = ['0' num2str(i-1)]; else, aux = num2str(i-1); end
% 
% 
%     Itex=imread([direc vis 'pose_' aux '.PNG']);
%     figure(1);
%     %title(['Reproyección de la pose ' num2str(i)])
%     subplot(121)
%     imshow(Itex)
%     hold on
%     plot(M(:,1,i,1),M(:,2,i,1),'xr')
%     plot(params.CameraParameters1.ReprojectedPoints(:,1,j),params.CameraParameters1.ReprojectedPoints(:,2,j),'+g')
%     title('Camera 1')
%     subplot(122)
%     imshow(img_ir{i},[])
%     hold on
%     plot(M(:,1,i,2),M(:,2,i,2),'xr')
%     plot(params.CameraParameters2.ReprojectedPoints(:,1,j),params.CameraParameters2.ReprojectedPoints(:,2,j),'+g')
%     title('Thermal camera')
%     hold off
%     sgtitle(['Pose ' aux])
%     pause(0.5); 
% 
% end
