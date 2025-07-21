n_poses=1;
load('D:\R3D\2023\multimodal\Experimento 21\Calib_TSC_0_3.mat'); % Directorio de archivo de calibración

for i=1:n_poses
    if i<11
        aux="0"+num2str(i-1);
    else
        aux=num2str(i-1);
    end
    load("D:\R3D\2023\multimodal\Experimento 23\camera 1\pose_"+aux+"\fp_phase.mat")
    Fy=phiy;
    Fx=phix;
    addpath('C:\Users\Eberto Benjumea\Desktop\OneDrive - Universidad Tecnológica de Bolívar\Doctorado\Proyectos\Codigos Reconstruccion 3D Learning\Codes Raul Purdue\StereoCalibCodes');


    suave=1;
    flag_rgb=0;
    Mask=ones(size(phiy));
    Mask=Mask(438:438+1200-1,584:584+1920-1);
    fx=Fx(438:438+1200-1,584:584+1920-1);
    fy=Fy(438:438+1200-1,584:584+1920-1);
    FringePitch=12;
    typeCalib=3;
    Fd='Fy';
    




    %Mask=ones(1024,1280);
    [XcM,YcM,ZcM,t1]=StereoReconstruction(fx,fy,Mask,FringePitch,StereoParams,typeCalib,Fd);
    ZcM=-ZcM;


    if suave==1
        %ZcM=imresize(ZcM,1,'bilinear');
        ZcM=remove_peaks(ZcM);
        %ZcM = imgaussfilt(ZcM);
        
        ZcM=medfilt2(ZcM,[9 9]);

    end
     figure;
%     masknan=isnan(ZcM);
%     imshow(masknan)
%     figure;
%     %M = median(A)
%     ZcM(masknan)=-500;
%     ZcM=medfilt2(ZcM,[11 11]);
    if flag_rgb==0
        s = surf(XcM,YcM,ZcM,'FaceColor', 'interp',...
                        'EdgeColor', 'none',...
                        'FaceLighting', 'phong');
    
        view(-160,80)
            set(gca, 'DataAspectRatio', [1, 1, 1])
            axis equal;
            view(0, 90);
            camlight right
            axis on
            grid on
            xlabel('x')
            ylabel('y')
            zlabel('z')
    else

        warp(XcM,YcM,ZcM,images);
    
        view(-160,80)
            set(gca, 'DataAspectRatio', [1, 1, 1])
            axis equal;
            view(0, 90);
            %camlight right
            %axis off
            grid off
            xlabel('x')
            ylabel('y')
            zlabel('z')
    end
    
    pause(1)
    if suave==0
        if i<10
            name_data=[Direc carp 'r3d_pose_0' num2str(i) '.mat'];
        else
            name_data=[Direc carp 'r3d_pose_' num2str(i) '.mat'];
        end
        save(name_data,'XcM','YcM','ZcM','images','t1')

    else
        if i<10
            name_data=[Direc carp 'r3d_smooth_pose_0' num2str(i) '.mat'];
        else
            name_data=[Direc carp 'r3d_smooth_pose_' num2str(i) '.mat'];
        end
        save(name_data,'XcM','YcM','ZcM','images','t1')

    end











end