close all
clc
clear all
direc='D:\R3D\Adquisiciones\Experimento 54\Digital features\';
figure;

correct_poses=[1 3:6 9:38 40:43 45:47 49:58 60:61 64:65 67:72 74:77 80];
for i=1:length(correct_poses)
    disp([num2str(i) ' of ' num2str(length(correct_poses))])
    if correct_poses(i)>9
        aux=num2str(correct_poses(i));
    else
        aux=['0' num2str(correct_poses(i))];
    end
    vis_texture1=imread([direc 'camera 1\pose_' aux '\im_51.png']);
    vis_texture2=imread([direc 'camera 2\pose_' aux '\im_51.png']);
    subplot(121)
    imshow(vis_texture1)
    subplot(122)
    imshow(vis_texture2)
    sgtitle(['Pose ' num2str(aux)])
    pause();
end