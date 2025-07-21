% Inputs:
% direc is the path where the main folder is stored.
% carp is the main folder name.
% NStep is the step number used in the creation of the fringe images.
% NBits is the number of bits.
% n_poses is the numer of poses
% m_imageWidth: image width
% m_imageHeight: image height
% FringePitch is set at 18. 

% Output:
% phase_camera is a cell with size n_poses x 2 where its content is
%  phase_camera=[{fx fy}; %For pose 1
%                {fx fy}; %For pose 2
%                {fx fy}; %For pose 3
%                 .   .  
%                 .   .  
%                 .   .  
%                {fx fy}]; %For pose n_poses
% By: OPI-lab Team
% More information: agmarrugo@utb.edu.co

function [phase_camera, MaskC]=estimate_phases(direc,carp,NStep,NBits,n_poses,m_imageWidth,m_imageHeight,format,inicio)

addpath('C:\Users\Eberto Benjumea\Desktop\OneDrive - Universidad Tecnológica de Bolívar\Doctorado\Proyectos\Codigos Reconstruccion 3D Learning\codigos calibración\code\Funciones_Hibrido')
%addpath('C:\Users\Beatriz\OneDrive - Universidad Tecnológica de Bolívar\Doctorado\Proyectos\Codigos Reconstruccion 3D Learning\codigos calibración\code\Funciones_Hibrido')

%carp = 'Brazo\';
%carp=camera1;
flag_color=0;
flag_roipoly = 0;
high_freq=0;    % 1:Yes and 0=No   
FringePitch = 18; %20
%NStep = 18; %10
%NBits = 7;
nitext = NStep+NBits+1;
phase_camera={};

%m_imageWidth = width; %1280; % 2048; %1280; 
%m_imageHeight = height; %1024; %1536; %1024;
%n_poses = 30;

Fx = zeros(m_imageHeight,m_imageWidth) ;
Fy = Fx; 
%images = Fx;
MaskC = Fx;
if flag_color==1
    images = zeros(m_imageWidth,m_imageHeight,3);
else
    images = Fx;
end

for i = inicio:n_poses
    Fx = zeros(m_imageHeight,m_imageWidth) ;
    Fy = Fx; 
    %images = Fx;
    MaskC = Fx;
    if(i<11), aux = ['0' num2str(i-1)]; else, aux = num2str(i-1); end

    % Deteccion de puntos - circulos
    Itex = imread([direc, carp ,'pose_',aux,'\H\P',num2str(FringePitch),'_F',num2str(nitext),format]);
    %images(:,:,i)= Itex;
    %images(:,:)= Itex(:,:,1);
    if flag_color==1
        ['D:\R3D\2023\Erik\Objetos\' 'color_pose_0' num2str(i-1) format]
        % Carga textura a color
        images= imread(['D:\R3D\2023\Erik\Objetos\' 'color_pose_0' num2str(i-1) format]);
    else
        images(:,:)= Itex(:,:,1);
    end
    %figure; imshow(images)
    %pause

    % Extraccion de fase
    if high_freq==0
        fileEx = [direc, carp ,'pose_',aux,'\V\P',num2str(FringePitch),'_F'];
        [Fx_aux, mask_h, texture_h, codeWord, ph] = UnwrapPhaseGrayCode_citedi(fileEx, FringePitch, NStep, NBits, [9 9],format);
        Fx(:,:) = Fx_aux;
        phx=ph;
    end

    fileEx = [direc, carp ,'pose_',aux,'\H\P',num2str(FringePitch),'_F'];
    [Fy_aux, mask_v, texture_v, codeWord, ph] = UnwrapPhaseGrayCode_citedi(fileEx, FringePitch, NStep, NBits, [9 9],format);
    Fy(:,:) = Fy_aux;
    phy=ph;
%     figure;
% 
%     imshow(ph)
    if(flag_roipoly)
        figure(81), imshow(Itex)
        title(['Imagen ',num2str(i)])
        rec  = double(roipoly);
        MaskC(:,:) = rec;
    else
        if high_freq==0
            MaskC(:,:) = mask_v.*mask_h;
        else
            MaskC(:,:) = mask_v;
        end
    end


    figure(80),imagesc(Fy_aux.*mask_v), title('Unwrapped H Phase - GrayCode');
    title(['Imagen ',num2str(i)])
    save([direc,carp,carp(1:end-1),'_pose_' aux '_phases.mat'],'Fx','Fy','images','MaskC')
    % For Erik
    %imwrite(Fy_aux.*mask_v,['for_mask_0' num2str(i-1) format])
    %imwrite(MaskC,[direc,carp, 'mask_pose_' aux format])
    %pause
    phase_camera{i,1}=Fx;
    phase_camera{i,2}=Fy;
end
end
