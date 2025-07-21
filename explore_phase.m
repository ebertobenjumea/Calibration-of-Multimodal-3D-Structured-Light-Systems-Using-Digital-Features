%     XcM=XcM(300:1000,300:1000);
%     YcM=YcM(300:1000,300:1000);
%     ZcM=ZcM(300:1000,300:1000);
    x=XcM(:);
    y=YcM(:);
    z=ZcM(:);
    x=x(~isnan(z));
        y=y(~isnan(z));
        z=z(~isnan(z));
    DM = [x, y, ones(size(z))];                             % Design Matrix
    %B = DM\z;                                               % Estimate Parameters
    B=lsqminnorm(DM,z);
    [X,Y] = meshgrid(linspace(min(x),max(x),50), linspace(min(y),max(y),50));
    Z = B(1)*x + B(2)*y + B(3)*ones(size(x));
    dif=z-Z;
    figure
    plot3(x,y,z)
title('z real')
    figure
    plot3(x,y,Z) 
    figure
    plot3(x,y,dif,'.b')
    zlim([-0.5 0.5])


RMSE = sqrt(mean((dif - 0).^2)); % Root Mean Squared Error

    figure;
    plot(isnan(x)-isnan(z))


    direct= [Direc 'Objects\pose_00\H\'];
    sum=double(zeros(1024,1280));
    figure; imshow(sum,[])
    for i=1:18
        if i<11
            img=imread([direct 'P18_F0' num2str(i-1) '.png']); % descomentar para archivos de raul
        else
            img=imread([direct 'P18_F' num2str(i-1) '.png']); % descomentar para archivos de raul
        end
        img=double(img);
        imshow(img,[])
        %img=imread([direct 'P18_F00']);
        sum=sum+img;
        pause(0.2);

    end
        figure; imshow(sum,[])
        title('sum')



        figure;
        texture=double(imread([direct 'P18_F26' '.png']));
        
                figure; imshow(texture,[])
                title('texture')



     figure;
     imagesc(Fx)
     figure; imagesc(Fy)
     fx_m=Fx.*Mask;
     fy_m=Fy.*Mask;
     figure;
     title('fx')
     imagesc(fx_m)
   figure;
     title('fy')
     imagesc(fy_m)

    %fx_m=fx_m(300:1000,300:1000);
    %fy_m=fy_m(300:1000,300:1000);

         [x,y]=meshgrid(1:size(fx_m,1),1:1024);
         x_nonan=x(:);
         y_nonan=y(:);
         fx_nonan=fx_m(:);
    fx_nonan(fx_nonan==0)=NaN;
    
    x_nonan=x_nonan(~isnan(fx_nonan));
    y_nonan=y_nonan(~isnan(fx_nonan));
    fx_nonan=fx_nonan(~isnan(fx_nonan));
    plot3(x_nonan,y_nonan,fx_nonan,'.b')

    DM = [x_nonan, y_nonan, ones(size(fx_nonan))];                             % Design Matrix
    %B = DM\z;                                               % Estimate Parameters
    B=lsqminnorm(DM,fx_nonan);
    %[X,Y] = meshgrid(linspace(min(x_nonan),max(x_nonan),50), linspace(min(y_nonan),max(y_nonan),50));
    Fx_nonan = B(1)*x_nonan + B(2)*y_nonan + B(3)*ones(size(x_nonan));
    dif=fx_nonan-Fx_nonan;
    figure;
    plot3(x_nonan,y_nonan,fx_nonan,'.b')
    title('Faase real')
    figure;
        plot3(x_nonan,y_nonan,Fx_nonan,'.b')
        title('Fase ideal')
            figure;
        plot3(x_nonan,y_nonan,dif,'.b')
        title('dif')


        %%
         [X,Y] = meshgrid(linspace(min(x_nonan),max(x_nonan),1000), linspace(min(y_nonan),max(y_nonan),1000));
        z_est_dig=griddata(x_nonan,y_nonan,dif,X,Y,'cubic');
figure;
s = surf(X,Y,z_est_dig,'FaceColor', 'interp',...
            'EdgeColor', 'none',...
            'FaceLighting', 'phong');


set(gca, 'DataAspectRatio', [1, 1, 1])
axis equal;
view(180,90);
camlight right
axis on
grid on
xlabel('x (mm)')
ylabel('y (mm)')
zlabel('z (mm)')

figure;
imagesc(z_est_dig)

figure;
plot(z_est_dig(25,:))





%%

    SE = strel("disk",5);

%Erosione la imagen con el elemento estructurante.

Mask= imerode(Mask,SE);

     figure;
     imagesc(Fx)
     figure; imagesc(Fy)
     fx_m=Fx.*Mask;
     fy_m=Fy.*Mask;
     figure;
     title('fx')
     imagesc(fx_m)
   figure;
     title('fy')
     imagesc(fy_m)



         [x,y]=meshgrid(1:1280,1:1024);
         x_nonan=x(:);
         y_nonan=y(:);
         fx_nonan=fy_m(:);
    fx_nonan(fx_nonan==0)=NaN;
    
    x_nonan=x_nonan(~isnan(fx_nonan));
    y_nonan=y_nonan(~isnan(fx_nonan));
    fx_nonan=fx_nonan(~isnan(fx_nonan));
    plot3(x_nonan,y_nonan,fx_nonan,'.b')

    DM = [x_nonan, y_nonan, ones(size(fx_nonan))];                             % Design Matrix
    %B = DM\z;                                               % Estimate Parameters
    B=lsqminnorm(DM,fx_nonan);
    %[X,Y] = meshgrid(linspace(min(x_nonan),max(x_nonan),50), linspace(min(y_nonan),max(y_nonan),50));
    Fx_nonan = B(1)*x_nonan + B(2)*y_nonan + B(3)*ones(size(x_nonan));
    dif=fx_nonan-Fx_nonan;
    figure;
    plot3(x_nonan,y_nonan,fx_nonan,'.b')
    title('Faase real')
    figure;
        plot3(x_nonan,y_nonan,Fx_nonan,'.b')
        title('Fase ideal')
            figure;
        plot3(x_nonan,y_nonan,dif,'.b')
        title('dif')


        %%
         [X,Y] = meshgrid(linspace(min(x_nonan),max(x_nonan),1000), linspace(min(y_nonan),max(y_nonan),1000));
        z_est_dig=griddata(x_nonan,y_nonan,dif,X,Y,'cubic');
figure;
s = surf(X,Y,z_est_dig,'FaceColor', 'interp',...
            'EdgeColor', 'none',...
            'FaceLighting', 'phong');


set(gca, 'DataAspectRatio', [1, 1, 1])
axis equal;
view(180,90);
camlight right
axis on
grid on
xlabel('x (mm)')
ylabel('y (mm)')
zlabel('z (mm)')

figure;
imagesc(z_est_dig)

figure;
plot(z_est_dig(25,:))
