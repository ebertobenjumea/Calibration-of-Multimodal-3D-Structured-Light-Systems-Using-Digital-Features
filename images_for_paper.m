%%
close all
clc
st=StereoParams;
cam1_params=st.CameraParameters1;
cam2_params=st.CameraParameters2;

%n_poses=3;
worldPoints={};
model_points={};
planes=[];
n=[];
p0_store=[];
H_stores={};

figure(1);
for i=1:n_poses
    disp(['Pose ' num2str(i-1)])
    mask_c2=~isnan(coor_c2{i});
    num_non_processed=sum(isnan(coor_c2{i}));
    disp(['Number of non processed points: ' num2str(num_non_processed(1))])
    coor_c2_aux=coor_c2{i}(mask_c2(:,1),:);
    coor_c1_aux=feature_points_c1_mirror{i}(mask_c2(:,1),:);
    coor_c1_undist = undistortPoints(double(coor_c1_aux),cam1_params);
    coor_c2_undist = undistortPoints(coor_c2_aux,cam2_params);
    %size(coor_c1_undist)
    %size(coor_c2_undist)
    worldPoints{i} = triangulate(coor_c1_undist,coor_c2_undist,StereoParams);
    XcM=worldPoints{i}(:,1);
    YcM=worldPoints{i}(:,2);
    ZcM=worldPoints{i}(:,3);
    

%    figure;
    subplot(131)
    plot3(XcM,YcM,ZcM,".")
    view([45 45])
    xlabel('x (mm)')
    ylabel('y (mm)')
    zlabel('z (mm)')
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
    
    subplot(132)
    Z = B(1)*x + B(2)*y + B(3)*ones(size(x));
    plot3(x,y,Z,'.r')
    view([45 45])
    %title('Reconstructed ideal 3D plane')
    xlabel('x (mm)')
    ylabel('y (mm)')
    zlabel('z (mm)')
    
    subplot(133)
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
    %title('Plane fitting error')
    

    pause();
%     close 
end


%%

close all
clc
st=StereoParams;
cam1_params=st.CameraParameters1;
cam2_params=st.CameraParameters2;

%n_poses=3;
worldPoints={};
model_points={};
planes=[];
n=[];
p0_store=[];
H_stores={};

figure(1);
for i=1:n_poses
    disp(['Pose ' num2str(i-1)])
    mask_c2=~isnan(coor_c2{i});
    num_non_processed=sum(isnan(coor_c2{i}));
    disp(['Number of non processed points: ' num2str(num_non_processed(1))])
    coor_c2_aux=coor_c2{i}(mask_c2(:,1),:);
    coor_c1_aux=feature_points_c1_mirror{i}(mask_c2(:,1),:);
    coor_c1_undist = undistortPoints(double(coor_c1_aux),cam1_params);
    coor_c2_undist = undistortPoints(coor_c2_aux,cam2_params);
    %size(coor_c1_undist)
    %size(coor_c2_undist)
    worldPoints{i} = triangulate(coor_c1_undist,coor_c2_undist,StereoParams);
    XcM=worldPoints{i}(:,1);
    YcM=worldPoints{i}(:,2);
    ZcM=worldPoints{i}(:,3);
    

%    figure;
    subplot(131)
    plot3(XcM,YcM,ZcM,".")
    %axis equal
    view([45 45])
    xlabel('x (mm)')
    ylabel('y (mm)')
    zlabel('z (mm)')
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
    
    subplot(132)
    Z = B(1)*x + B(2)*y + B(3)*ones(size(x));
    plot3(x,y,Z,'.r')
    %axis equal
    view([45 45])
    %title('Reconstructed ideal 3D plane')
    xlabel('x (mm)')
    ylabel('y (mm)')
    zlabel('z (mm)')
    
    %subplot(132)
    error_z=z-Z;
    %plot3(x,y,error_z,'.r')
%     dif_p=[x-x y-y z-Z];
%     error=[];
%     for j=1:size(dif_p,1)
%         error(j)=norm(dif_p(j,:));
%     end
%     plot3(x,y,error,'.r')
    %view([45 45])
    %zlim([-2 2])
    %xlabel('x (mm)')
    %ylabel('y (mm)')
   % zlabel('z (mm)')
    %title('Plane fitting error')
    

    %subplot(234)
    x0=x(1);
    y0=y(1);
    z0=Z(1);
    norm_plan=n(:,i);
    nz=norm_plan./norm(norm_plan);
    %figure; 
    ps=search_x_axis(x,y,Z);
    p_0=ps(:,2);%+[10;-10];
    p_1=ps(:,1);
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
    %xlabel('x (mm)')
    %ylabel('y (mm)')
    %zlabel('z (mm)')
    %title('x and y coordinates')
    
    subplot(133)
    plot3(model_points{i}(1,:),model_points{i}(2,:),model_points{i}(3,:),'.g')
    %axis equal
    zlim([-1 1])
    xlabel('x (mm)')
    ylabel('y (mm)')
    zlabel('z (mm)')
        view([45 45])
    %title('Points on model coordinate system')
    %sgtitle(['Pose ' num2str(i-1)])
    pause();
%     close 
end