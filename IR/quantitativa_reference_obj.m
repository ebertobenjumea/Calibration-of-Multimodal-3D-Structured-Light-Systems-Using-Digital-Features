close all
clc
clear all
direc='D:\R3D\Adquisiciones\Experimento 54\Objects\SL\';
aux=num2str(6);
% Cargamos digital features image
pc_dig = pcread([direc 'Digital Features 1 UTB\obj_vis_ir_0' aux '.ply']);
figure;
subplot(121)
pcshow(pc_dig);
xlabel('x (mm)')
ylabel('y (mm)')
zlabel('z (mm)')
c = colorbar; 
c.Label.String = '° C';
set(gca, 'Zdir', 'reverse')
set(gca, 'Xdir', 'reverse')
view(30,45)
grid off
zlim([350 430])
title('Digital features I')

% Cargamos traditional image
pc_trad = pcread([direc 'Two Stages 2\obj_vis_ir_0' aux '.ply']);
subplot(122)
pcshow(pc_trad);
xlabel('x (mm)')
ylabel('y (mm)')
zlabel('z (mm)')
c = colorbar; 
c.Label.String = '° C';
set(gca, 'Zdir', 'reverse')
set(gca, 'Xdir', 'reverse')
view(30,45)

grid off
zlim([350 430])
title('Traditional')



%
%% Creamos matrices de información
% Digital features I
x_dig=pc_dig.Location(:,1);
y_dig=pc_dig.Location(:,2);
z_dig=pc_dig.Location(:,3);
t_dig=pc_dig.Intensity;
%[X_dig,Y_dig]=meshgrid([min(x_dig):0.05:max(x_dig)], [min(y_dig):0.05:max(y_dig)]);

[X_dig,Y_dig]=meshgrid(linspace(min(x_dig),max(x_dig),1280), linspace(min(y_dig),max(y_dig),1024));



z_est_dig=griddata(x_dig,y_dig,z_dig,X_dig,Y_dig,'cubic');
t_est_dig=griddata(x_dig,y_dig,t_dig,X_dig,Y_dig,'cubic');
z_est_dig(z_est_dig>max(z_dig))=NaN;
t_est_dig(t_est_dig>max(t_dig))=NaN;
z_est_dig(z_est_dig<min(z_dig))=NaN;
t_est_dig(t_est_dig<min(t_dig))=NaN;

% Visualización de imagenes 3D y térmicas
figure;
subplot(121)
imagesc(z_est_dig)
colorbar
axis equal
subplot(122)
imagesc(t_est_dig)
colorbar
axis equal

% Visualizacion 3D
figure;
s = surf(-X_dig,Y_dig,-z_est_dig,'FaceColor', 'interp',...
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




disp('Cantidad de datos NaN:')
sum(sum(isnan(z_est_dig))) % Ojo con esto
sum(sum(isnan(t_est_dig))) % Ojo con esto

% aux_t_est=t_est;
% mean_t_est=mean(mean(t_est,"omitnan"),"omitnan");
% aux_t_est(isnan(t_est))=mean_t_est;
% 
% aux_z_est=z_est;
% mean_z_est=mean(mean(z_est,"omitnan"),"omitnan");
% aux_z_est(isnan(z_est))=mean_z_est;
% 
% sum(sum(isnan(aux_z_est))) % Ojo con esto
% sum(sum(isnan(aux_t_est))) % Ojo con esto
% save('datos_mapeo.mat',"aux_t_est","aux_z_est","X","Y")

% Traditional
x_trad=pc_trad.Location(:,1);
y_trad=pc_trad.Location(:,2);
z_trad=pc_trad.Location(:,3);
t_trad=pc_trad.Intensity;
%[X_trad,Y_trad]=meshgrid([min(x_trad):0.05:max(x_trad)], [min(y_trad):0.05:max(y_trad)]);
[X_trad,Y_trad]=meshgrid(linspace(min(x_trad),max(x_trad),1280), linspace(min(y_trad),max(y_trad),1024));



z_est_trad=griddata(x_trad,y_trad,z_trad,X_trad,Y_trad,'cubic');
t_est_trad=griddata(x_trad,y_trad,t_trad,X_trad,Y_trad,'cubic');
z_est_trad(z_est_trad>max(z_trad))=NaN;
t_est_trad(t_est_trad>max(t_trad))=NaN;
z_est_trad(z_est_trad<min(z_trad))=NaN;
t_est_trad(t_est_trad<min(t_trad))=NaN;

% Visualización de imagenes 3D y térmicas
figure;
subplot(121)
imagesc(z_est_trad)
colorbar
axis equal
subplot(122)
imagesc(t_est_trad)
colorbar
axis equal

% Visualizacion 3D
figure;
s = surf(-X_trad,Y_trad,-z_est_trad,'FaceColor', 'interp',...
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


disp('Cantidad de datos NaN:')
sum(sum(isnan(z_est_trad))) % Ojo con esto
sum(sum(isnan(t_est_trad))) % Ojo con esto


%%


figure;
s = surf(-X_trad(100:end-100,100:end-100),Y_trad(100:end-100,100:end-100),-z_est_trad(100:end-100,100:end-100),'FaceColor', 'interp',...
            'EdgeColor', 'none',...
            'FaceLighting', 'phong');


set(gca, 'DataAspectRatio', [1, 1, 1])
axis equal;
view(180,90);
%camlight right
axis on
grid on
xlabel('x (mm)')
ylabel('y (mm)')
zlabel('z (mm)')
%figure;
hold on

s = surf(-X_dig(100:end-100,100:end-100),Y_dig(100:end-100,100:end-100),-z_est_dig(100:end-100,100:end-100),'FaceColor', 'b',...
            'EdgeColor', 'none',...
            'FaceLighting', 'phong');
set(gca, 'DataAspectRatio', [1, 1, 1])

%%
close all

figure;

s = surf(-X_trad,Y_trad,-z_est_trad,'FaceColor', 'interp',...
            'EdgeColor', 'none',...
            'FaceLighting', 'phong');


set(gca, 'DataAspectRatio', [1, 1, 1])
axis equal;
view(180,90);
%camlight right
axis on
grid on
xlabel('x (mm)')
ylabel('y (mm)')
zlabel('z (mm)')
%figure;
hold on
title('Trad vs Dig 2')



s = surf(-X_dig2,Y_dig2,-z_est_dig2,'FaceColor', 'r',...
            'EdgeColor', 'none',...
            'FaceLighting', 'phong');
set(gca, 'DataAspectRatio', [1, 1, 1])


figure;

s = surf(-X_dig,Y_dig,-z_est_dig,'FaceColor', 'b',...
            'EdgeColor', 'none',...
            'FaceLighting', 'phong');
set(gca, 'DataAspectRatio', [1, 1, 1])
hold on

s = surf(-X_dig2,Y_dig2,-z_est_dig2,'FaceColor', 'r',...
            'EdgeColor', 'none',...
            'FaceLighting', 'phong');
set(gca, 'DataAspectRatio', [1, 1, 1])

axis equal;
view(180,90);
%camlight right
axis on
grid on
xlabel('x (mm)')
ylabel('y (mm)')
zlabel('z (mm)')
%figure;
hold on
title('Dig vs Dig 2')

figure;
s = surf(-X_trad,Y_trad,-z_est_trad,'FaceColor', 'interp',...
            'EdgeColor', 'none',...
            'FaceLighting', 'phong');


set(gca, 'DataAspectRatio', [1, 1, 1])
axis equal;
view(180,90);
%camlight right
axis on
grid on
xlabel('x (mm)')
ylabel('y (mm)')
zlabel('z (mm)')
%figure;
hold on

s = surf(-X_dig,Y_dig,-z_est_dig,'FaceColor', 'b',...
            'EdgeColor', 'none',...
            'FaceLighting', 'phong');
set(gca, 'DataAspectRatio', [1, 1, 1])
hold on

s = surf(-X_dig2,Y_dig2,-z_est_dig2,'FaceColor', 'r',...
            'EdgeColor', 'none',...
            'FaceLighting', 'phong');
set(gca, 'DataAspectRatio', [1, 1, 1])
%% Vista de perfiles horizontales - Digital Features
% Digital Features
close all
n_profile=300;

figure;
imagesc(z_est_dig,[360 430])
yline(n_profile,'--k')
colorbar

figure;
imagesc(t_est_dig)
yline(n_profile,'--k')
colorbar



tam=length(X_dig(n_profile,:));
largo_profile=120:tam-120;
figure(6);
subplot(211)
plot(z_est_dig(n_profile,largo_profile))
xlabel('X (mm)')
ylabel('Z (mm)')
hold on
subplot(212)
plot(t_est_dig(n_profile,largo_profile))
xlabel('X (mm)')
ylabel('t (°C)')
hold on

mean_dig=mean(z_est_dig(n_profile,largo_profile),'omitnan');
aux_dig=z_est_dig(n_profile,largo_profile)-mean_dig;
% Vista de perfiles en la misma grafica sin perfil primario
profile=-aux_dig;
profile=profile(~isnan(profile));
profile_z_dig=profile;
X_dig_sin_nan=X_dig(n_profile,:);
X_dig_sin_nan=X_dig_sin_nan(~isnan(profile));


p = polyfit(1:length(profile),profile,1); 
f = polyval(p,1:length(profile)); 
profile_z_dig_sin_primary=profile-f;
% figure;
% plot(f)
% hold on
% plot(profile)

%quitamos perfil primario
%Buscamos el perfil termico
profile_t_dig=t_est_dig(n_profile,largo_profile);
profile_t_dig=profile_t_dig(~isnan(profile));

size(profile_z_dig_sin_primary)
size(profile_t_dig)
figure(7);

yyaxis left
offset=-min(profile_z_dig_sin_primary);
offset=0;
plot((profile_z_dig_sin_primary+1*offset),'b')
ylabel('Z (mm)')
% yt = get(gca, 'YTick');
% set(gca, 'YTick',yt, 'ZTickLabel',fliplr(yt))
%ylim([-1 1])
yyaxis right
plot(profile_t_dig)
xlabel('X (mm)')
ylabel('t (°C)')
title('Profile in x')
hold on


%quitamos perfil secundario
p1 = polyfit(1:length(profile_z_dig_sin_primary),profile_z_dig_sin_primary,2); 
f1 = polyval(p1,1:length(profile_z_dig_sin_primary)); 
figure(8);
subplot(311)
yyaxis left
plot((profile_z_dig_sin_primary-f1))
ylabel('Z (mm)')
%ylim([-1 1])
yyaxis right
plot(profile_t_dig)
xlabel('X (mm)')
ylabel('t (°C)')
title('Digital')
grid on


%figure zoom
figure(9);

aux=(profile_z_dig_sin_primary-f1);
profile_final_z_dig=aux;
tam_aux=1:length(aux);
%tam_aux=550:700;
%tam_aux=570:590


yyaxis left
offset=-min(aux(tam_aux));
plot(aux(tam_aux)+offset)
ylabel('Z (mm)')
%ylim([-1 1])
yyaxis right
plot(profile_t_dig(tam_aux))
xlabel('X (px)')
ylabel('t (°C)')
%title('Digital')
grid on
hold on


%% Vista de perfiles horizontales - Traditional
% Digital Features
%close all
n_profile=300;

figure;
imagesc(z_est_trad,[360 430])
yline(n_profile,'--k')
colorbar

figure;
imagesc(t_est_trad)
yline(n_profile,'--k')
colorbar


tam=length(X_trad(n_profile,:));
largo_profile=120:tam-120;
y_linea=n_profile;
figure(6);
subplot(211)
plot(z_est_trad(n_profile,largo_profile))
xlabel('X (mm)')
ylabel('Z (mm)')
legend
subplot(212)
plot(t_est_trad(n_profile,largo_profile))
xlabel('X (mm)')
ylabel('t (°C)')
legend


mean_trad=mean(z_est_trad(n_profile,largo_profile),'omitnan');
aux_trad=z_est_trad(n_profile,largo_profile)-mean_trad;
% Vista de perfiles en la misma grafica sin perfil primario
profile=-aux_trad;

% Vista de perfiles en la misma grafica sin perfil primario
profile=-z_est_trad(n_profile,largo_profile);
profile=profile(~isnan(profile));
profile_z_trad=profile;
X_trad_sin_nan=X_trad(n_profile,:);
X_trad_sin_nan=X_trad_sin_nan(~isnan(profile));


p = polyfit(1:length(profile),profile,1); 
f = polyval(p,1:length(profile)); 
profile_z_trad_sin_primary=profile-f;
% figure;
% plot(f)
% hold on
% plot(profile)

%quitamos perfil primario
%Buscamos el perfil termico
profile_t_trad=t_est_trad(n_profile,largo_profile);
profile_t_trad=profile_t_trad(~isnan(profile));

size(profile_z_trad_sin_primary)
size(profile_t_trad)
figure(7);

yyaxis left
offset=-min(profile_z_trad_sin_primary);
offset=0;
plot((profile_z_trad_sin_primary+1*offset),'b')
ylabel('Z (mm)')
% yt = get(gca, 'YTick');
% set(gca, 'YTick',yt, 'ZTickLabel',fliplr(yt))
%ylim([-1 1])
yyaxis right
plot(profile_t_trad)
xlabel('X (mm)')
ylabel('t (°C)')
title('Profile in x')



%quitamos perfil secundario
p1 = polyfit(1:length(profile_z_trad_sin_primary),profile_z_trad_sin_primary,2); 
f2 = polyval(p1,1:length(profile_z_trad_sin_primary)); 
figure(8);
subplot(212)
yyaxis left
plot((profile_z_trad_sin_primary-f2))
ylabel('Z (mm)')
%ylim([-1 1])
yyaxis right
plot(profile_t_trad)
xlabel('X (mm)')
ylabel('t (°C)')
title('Traditional')
grid on




figure(9);
aux=(profile_z_trad_sin_primary-f2);
profile_final_z_trad=aux;
offset=-min(profile_final_z_trad);
yyaxis left
plot(aux(tam_aux)+offset,'-k')
ylabel('Z (mm)')
%ylim([-1 1])
yyaxis right
plot(profile_t_trad(tam_aux),'-g')
xlabel('X (px)')
ylabel('t (°C)')
%title('Digital')
grid on
legend('Z  by Digital features 1','Z  by Traditional','t by Digital features 1','t by Traditional')
grid minor
%%
mean_z_dig=mean(profile_z_dig)
RMSE_z_dig=sqrt(sum(((profile_z_dig-mean_z_dig).^2)))/(length(profile_z_dig))
mean_z_trad=mean(profile_z_trad)
RMSE_z_trad=sqrt(sum(((profile_z_trad-mean_z_trad).^2)))/(length(profile_z_trad))

mean_t_dig=mean(profile_t_dig)
RMSE_t_dig=sqrt(sum(((profile_t_dig-mean_t_dig).^2)))/(length(profile_t_dig))
mean_t_trad=mean(profile_t_trad)
RMSE_t_trad=sqrt(sum(((profile_t_trad-mean_t_trad).^2)))/(length(profile_t_trad))


%% Entre señales
mean_z_dig=mean(profile_z_dig)
%RMSE_z_dig=sqrt(sum(((profile_z_dig-mean_z_dig).^2)))/(length(profile_z_dig))
mean_z_trad=mean(profile_z_trad)
RMSE_z=sqrt(sum(((profile_z_trad-mean_z_trad-(profile_z_dig-mean_z_dig)).^2)))/(length(profile_z_trad))

mean_t_dig=mean(profile_t_dig)
mean_t_trad=mean(profile_t_trad)
RMSE_t_=sqrt(sum(((profile_t_trad-mean_t_trad-(profile_t_dig-mean_t_dig)).^2)))/(length(profile_t_trad))


%% Entre imagenes
figure;
subplot(221)
imagesc(z_est_trad)
title('Z trad')
subplot(222)
imagesc(z_est_dig)
title('Z dig')
subplot(223)
imagesc(t_est_trad)
title('t trad')
subplot(224)
imagesc(t_est_dig)
title('t dig')

%%
close all
t_est_trad_copy=t_est_trad;
figure; 
imshow(t_est_trad,[])
B = sort(t_est_trad(:));
C = unique(B);
t_est_trad(isnan(t_est_trad))=C(2); %mean(mean(t_sin_perfil,"omitnan"),"omitnan");
t_est_trad(t_est_trad==0)=C(2); %mean(mean(t_sin_perfil,"omitnan"),"omitnan");

X=t_est_trad;
Xmin = min(X(:));
Xmax = max(X(:));
if isequal(Xmax,Xmin)
    X = 0*X;
else
    X = (X - Xmin) ./ (Xmax - Xmin);
end

% Threshold image - adaptive threshold
BW = imbinarize(X, 'adaptive', 'Sensitivity', 0.500000, 'ForegroundPolarity', 'bright');
figure; imshow(BW)
% Open mask with disk
radius = 10;
decomposition = 0;
se = strel('disk', radius, decomposition);
a = imopen(BW, se);

% Threshold image - adaptive threshold
figure; imshow(a)
[mask_w_w,numWhite] = bwlabel(a);
mask_w_w = bwareaopen(mask_w_w,10000);  %Retirando manchas blancas
figure; imshow(mask_w_w)

%a=(imbinarize(t_sin_perfil));
a=mask_w_w ;
b= ~bwareaopen(a,25000);   %Retirando manchas negras
a=logical(a.*b);
 figure; 
 imshow(a,[])
% [puntos, radii, metric] = imfindcircles(a,[50 70],'ObjectPolarity','bright',Sensitivity=0.9);
% figure;
% imagesc(a)
% hold on
% plot(puntos(:,1),puntos(:,2),'ok')
% viscircles(puntos,radii)

figure;
imagesc(a)
stats = regionprops("table",a,"Centroid", ...
    "MajorAxisLength","MinorAxisLength")
centers_t = stats.Centroid;
diameters = mean([stats.MajorAxisLength stats.MinorAxisLength],2);
radii_t = diameters/2;
hold on
viscircles(centers_t,radii_t)
plot(centers_t(:,1),centers_t(:,2),'+k')

%%
close all
t_est_dig_copy=t_est_dig;
figure; 
imshow(t_est_dig,[])
B = sort(t_est_dig(:));
C = unique(B);
t_est_dig(isnan(t_est_dig))=C(2); %mean(mean(t_sin_perfil,"omitnan"),"omitnan");
t_est_dig(t_est_dig==0)=C(2); %mean(mean(t_sin_perfil,"omitnan"),"omitnan");

X=t_est_dig;
Xmin = min(X(:));
Xmax = max(X(:));
if isequal(Xmax,Xmin)
    X = 0*X;
else
    X = (X - Xmin) ./ (Xmax - Xmin);
end

% Threshold image - adaptive threshold
BW = imbinarize(X, 'adaptive', 'Sensitivity', 0.500000, 'ForegroundPolarity', 'bright');
figure; imshow(BW)
% Open mask with disk
radius = 10;
decomposition = 0;
se = strel('disk', radius, decomposition);
a = imopen(BW, se);

% Threshold image - adaptive threshold
figure; imshow(a)
[mask_w_w,numWhite] = bwlabel(a);
mask_w_w = bwareaopen(mask_w_w,10000);  %Retirando manchas blancas
figure; imshow(mask_w_w)

%a=(imbinarize(t_sin_perfil));
a=mask_w_w ;
b= ~bwareaopen(a,12000);   %Retirando manchas negras
a=logical(a.*b);
 figure; 
 imshow(a,[])
% [puntos, radii, metric] = imfindcircles(a,[50 70],'ObjectPolarity','bright',Sensitivity=0.9);
% figure;
% imagesc(a)
% hold on
% plot(puntos(:,1),puntos(:,2),'ok')
% viscircles(puntos,radii)

figure;
imagesc(a)
stats = regionprops("table",a,"Centroid", ...
    "MajorAxisLength","MinorAxisLength")
centers_dig = stats.Centroid;
diameters = mean([stats.MajorAxisLength stats.MinorAxisLength],2);
radii_dig = diameters/2;
hold on
viscircles(centers_dig,radii_dig)
plot(centers_dig(:,1),centers_dig(:,2),'+k')



close all
figure;
imagesc(t_est_trad_copy)
hold on
for i=1:size(centers_t,1)
    plot(centers_t(i,1),centers_t(i,2),'+b')
    plot(centers_dig(i,1),centers_dig(i,2),'xk')
    %rms_t=norm(centers_t-centers_dig)
    pause(0.1);
end
rms_t=(sum(vecnorm(centers_t-centers_dig,2,2)))/size(centers_t,1);


errormapt=NaN.*ones(size(t_est_trad_copy));
error_pixel=vecnorm(centers_t-centers_dig,2,2);
for i=1:size(centers_t,1)
    errormapt(floor(centers_t(i,2)),floor(centers_t(i,1)))=error_pixel(i);
end
figure; imagesc(errormapt)
colorbar

figure;
stem(error_pixel,'ob')


[x_error,y_error]=meshgrid(1:1280,1:1024);
x_error=x_error(:);
y_error=y_error(:);
errormapt_interpol= griddata(double(centers_t(:,1)),double(centers_t(:,2)),error_pixel,double(x_error),double(y_error),"cubic");
figure; imagesc(reshape(errormapt_interpol,[1024, 1280]))
colorbar


figure;
imagesc(t_est_trad_copy-t_est_dig_copy)
colorbar
%% DIG 2


%%

t_est_dig_copy=t_est_dig2;
figure; 
imshow(t_est_dig2,[])
B = sort(t_est_dig2(:));
C = unique(B);
t_est_dig2(isnan(t_est_dig2))=C(2); %mean(mean(t_sin_perfil,"omitnan"),"omitnan");
t_est_dig2(t_est_dig2==0)=C(2); %mean(mean(t_sin_perfil,"omitnan"),"omitnan");

X=t_est_dig2;
Xmin = min(X(:));
Xmax = max(X(:));
if isequal(Xmax,Xmin)
    X = 0*X;
else
    X = (X - Xmin) ./ (Xmax - Xmin);
end

% Threshold image - adaptive threshold
BW = imbinarize(X, 'adaptive', 'Sensitivity', 0.500000, 'ForegroundPolarity', 'bright');
figure; imshow(BW)
% Open mask with disk
radius = 10;
decomposition = 0;
se = strel('disk', radius, decomposition);
a = imopen(BW, se);

% Threshold image - adaptive threshold
figure; imshow(a)
[mask_w_w,numWhite] = bwlabel(a);
mask_w_w = bwareaopen(mask_w_w,10000);  %Retirando manchas blancas
figure; imshow(mask_w_w)

%a=(imbinarize(t_sin_perfil));
a=mask_w_w ;
b= ~bwareaopen(a,12000);   %Retirando manchas negras
a=logical(a.*b);
 figure; 
 imshow(a,[])
% [puntos, radii, metric] = imfindcircles(a,[50 70],'ObjectPolarity','bright',Sensitivity=0.9);
% figure;
% imagesc(a)
% hold on
% plot(puntos(:,1),puntos(:,2),'ok')
% viscircles(puntos,radii)

figure;
imagesc(a)
stats = regionprops("table",a,"Centroid", ...
    "MajorAxisLength","MinorAxisLength")
centers_dig2 = stats.Centroid;
diameters = mean([stats.MajorAxisLength stats.MinorAxisLength],2);
radii_dig = diameters/2;
hold on
viscircles(centers_dig2,radii_dig)
plot(centers_dig2(:,1),centers_dig2(:,2),'+k')




figure;
imagesc(t_est_trad_copy)
hold on
for i=1:size(centers_t,1)
    plot(centers_t(i,1),centers_t(i,2),'+b')
    plot(centers_dig2(i,1),centers_dig2(i,2),'xk')
    %rms_t=norm(centers_t-centers_dig)
    pause(0.1);
end
rms_t=(sum(vecnorm(centers_t-centers_dig2,2,2)))/size(centers_t,1);


errormapt=NaN.*ones(size(t_est_trad_copy));
error_pixel=vecnorm(centers_t-centers_dig2,2,2);
for i=1:size(centers_t,1)
    errormapt(floor(centers_t(i,2)),floor(centers_t(i,1)))=error_pixel(i);
end
figure; imagesc(errormapt)
colorbar

figure;
stem(error_pixel,'ob')


[x_error,y_error]=meshgrid(1:1280,1:1024);
x_error=x_error(:);
y_error=y_error(:);
errormapt_interpol= griddata(double(centers_t(:,1)),double(centers_t(:,2)),error_pixel,double(x_error),double(y_error),"cubic");
figure; imagesc(reshape(errormapt_interpol,[1024, 1280]))
colorbar


figure;
imagesc(t_est_trad_copy-t_est_dig_copy)
colorbar
%%
close all
figure;
%yyaxis left
plot(-z_est_dig(n_profile,largo_profile),'b')
hold on
plot(-z_est_trad(n_profile,largo_profile),'r')
xlabel('Pixels')
ylabel('Z (mm)')
legend({'Digital Features','Traditional'})
grid on


figure;

plot(t_est_dig(n_profile,largo_profile),'g')
hold on
plot(t_est_trad(n_profile,largo_profile),'k')
xlabel('Pixels')
legend({'Digital Features','Traditional'})
ylabel('t (°C)')
grid on

figure;
yyaxis left
plot((profile_z_dig_sin_primary-f1),'b')
ylabel('Z (mm)')
yyaxis right
plot(t_est_dig(n_profile,largo_profile),'g')
hold on
plot(t_est_trad(n_profile,largo_profile),'k')
ylabel('t (°C)')
xlabel('Pixels')
legend({'Z','t by Digital features','t by Traditional'})


%%
close all

xx=1:1000;
yy=(-50/1279)*xx+735;
figure;
imagesc(t_est_dig)
hold on
plot(xx,yy,'.k')

figure;
imagesc(z_est_dig,[-450 -395])
colorbar
hold on
plot(xx,yy,'.k')

z_est_diag=[];
t_est_diag=[];
for i=1:1280
    z_est_diag(i)=-z_est_dig(round(yy(i)),xx(i));
    t_est_diag(i)=t_est_dig(round(yy(i)),xx(i));
end

%z_est_diag=z_est_dig(xx,yy);
%t_est_diag=t_est_dig(xx,yy);
figure;
plot(z_est_diag)
figure;
plot(t_est_diag)





tam=length(X_dig(n_profile,:));
largo_profile=350:tam-120;
figure;
subplot(211)
plot(z_est_diag(largo_profile))
xlabel('X (mm)')
ylabel('Z (mm)')
subplot(212)
plot(t_est_diag(largo_profile))
xlabel('X (mm)')
ylabel('t (°C)')


% Vista de perfiles en la misma grafica sin perfil primario
profile=z_est_diag(largo_profile);
profile=profile(~isnan(profile));
X_dig_sin_nan=X_dig(n_profile,:);
X_dig_sin_nan=X_dig_sin_nan(~isnan(profile));


p = polyfit(1:length(profile),profile,1); 
f = polyval(p,1:length(profile)); 
profile_z_diag_sin_primary=profile-f;
% figure;
% plot(f)
% hold on
% plot(profile)

%quitamos perfil primario
%Buscamos el perfil termico
profile_t_diag=t_est_diag(largo_profile);
profile_t_diag=profile_t_diag(~isnan(profile));

size(profile_z_diag_sin_primary)
size(profile_t_diag)
figure(7);
subplot(211)
yyaxis left
offset=-min(profile_z_diag_sin_primary);
%offset=0;
plot((profile_z_diag_sin_primary+1*offset),'b')
ylabel('Z (mm)')
% yt = get(gca, 'YTick');
% set(gca, 'YTick',yt, 'ZTickLabel',fliplr(yt))
%ylim([-1 1])
yyaxis right
plot(profile_t_diag)
xlabel('X (mm)')
ylabel('t (°C)')
title('Profile in x')



%quitamos perfil secundario
p1 = polyfit(1:length(profile_z_diag_sin_primary),profile_z_diag_sin_primary,3); 
f1 = polyval(p1,1:length(profile_z_diag_sin_primary)); 
figure(8);
subplot(211)
yyaxis left
plot((profile_z_diag_sin_primary-f1))
ylabel('Z (mm)')
%ylim([-1 1])
yyaxis right
plot(profile_t_diag)
xlabel('X (mm)')
ylabel('t (°C)')
title('Digital')
grid on


%%

close all
figure;
%yyaxis left
plot(-z_est_dig(n_profile,largo_profile),'b')
hold on
plot(-z_est_trad(n_profile,largo_profile),'r')
xlabel('Pixels')
ylabel('Z (mm)')
legend({'Digital Features','Traditional'})
grid on


figure;

plot(t_est_dig(n_profile,largo_profile),'g')
hold on
plot(t_est_trad(n_profile,largo_profile),'k')
xlabel('Pixels')
legend({'Digital Features','Traditional'})
ylabel('t (°C)')
grid on

figure;
yyaxis left
plot((profile_z_dig_sin_primary-f1),'b')
ylabel('Z (mm)')
yyaxis right
plot(t_est_dig(n_profile,largo_profile),'g')
hold on
plot(t_est_trad(n_profile,largo_profile),'k')
ylabel('t (°C)')
xlabel('Pixels')
legend({'Z','t by Digital features','t by Traditional'})

%% Tradicional
xx=350:1000;
yy=(-50/1279)*xx+735;
figure;
imagesc(t_est_trad)
hold on
plot(xx,yy,'.k')

z_est_diag=[];
t_est_diag=[];
x_plt_trad=[];
for i=1:length(xx)
    z_est_diag(i)=z_est_trad(round(yy(i)),xx(i));
    t_est_diag(i)=t_est_trad(round(yy(i)),xx(i));
    x_plt_trad(i)=X_trad(round(yy(i)),xx(i));
end

%z_est_diag=z_est_dig(xx,yy);
%t_est_diag=t_est_dig(xx,yy);




figure;
subplot(211)
plot(z_est_diag)
xlabel('X (mm)')
ylabel('Z (mm)')
subplot(212)
plot(t_est_diag)
xlabel('X (mm)')
ylabel('t (°C)')

% Vista de perfiles en la misma grafica sin perfil primario
profile=z_est_diag;
profile=profile(~isnan(profile));
X_dig_sin_nan=xx;
X_dig_sin_nan=X_dig_sin_nan(~isnan(profile));


p = polyfit(1:length(profile),profile,1); 
f = polyval(p,1:length(profile)); 
profile_z_diag_sin_primary=profile-f;
% figure;
% plot(f)
% hold on
% plot(profile)

%quitamos perfil primario
%Buscamos el perfil termico
profile_t_diag=t_est_diag;
profile_t_diag=profile_t_diag(~isnan(profile));

size(profile_z_diag_sin_primary)
size(profile_t_diag)
figure(7);
subplot(212)
yyaxis left
offset=-min(profile_z_diag_sin_primary);
%offset=0;
plot((profile_z_diag_sin_primary+1*offset),'b')
ylabel('Z (mm)')
% yt = get(gca, 'YTick');
% set(gca, 'YTick',yt, 'ZTickLabel',fliplr(yt))
%ylim([-1 1])
yyaxis right
plot(profile_t_diag)
xlabel('X (mm)')
ylabel('t (°C)')
title('Profile in x')



%quitamos perfil secundario
p1 = polyfit(1:length(profile_z_diag_sin_primary),profile_z_diag_sin_primary,3); 
f1 = polyval(p1,1:length(profile_z_diag_sin_primary)); 
figure(8);
subplot(212)
yyaxis left
plot(x_plt_trad,-(profile_z_diag_sin_primary-f1))
ylabel('Z (mm)')
%ylim([-1 1])
yyaxis right
plot(x_plt_trad,profile_t_diag)
xlabel('Profile (mm)')
ylabel('t (°C)')
%title('Trad')
grid on
%% Digital

xx=350:1000;
yy=(-50/1279)*xx+735;
figure;
imagesc(t_est_dig)
hold on
plot(xx,yy,'.k')

z_est_diag=[];
t_est_diag=[];
x_plt_dig=[];
for i=1:length(xx)
    z_est_diag(i)=z_est_dig(round(yy(i)),xx(i));
    t_est_diag(i)=t_est_dig(round(yy(i)),xx(i));
        x_plt_dig(i)=X_dig(round(yy(i)),xx(i));
end

%z_est_diag=z_est_dig(xx,yy);
%t_est_diag=t_est_dig(xx,yy);




figure;
subplot(211)
plot(z_est_diag)
xlabel('X (mm)')
ylabel('Z (mm)')
subplot(212)
plot(t_est_diag)
xlabel('X (mm)')
ylabel('t (°C)')

% Vista de perfiles en la misma grafica sin perfil primario
profile=z_est_diag;
profile=profile(~isnan(profile));
X_dig_sin_nan=xx;
X_dig_sin_nan=X_dig_sin_nan(~isnan(profile));


p = polyfit(1:length(profile),profile,1); 
f = polyval(p,1:length(profile)); 
profile_z_diag_sin_primary=profile-f;

figure;
plot(profile+400)
hold on
plot(f+400)
hold on
plot(profile_z_diag_sin_primary)
% figure;
% plot(f)
% hold on
% plot(profile)

%quitamos perfil primario
%Buscamos el perfil termico
profile_t_diag=t_est_diag;
profile_t_diag=profile_t_diag(~isnan(profile));

size(profile_z_diag_sin_primary)
size(profile_t_diag)
figure(7);
subplot(211)
yyaxis left
offset=-min(profile_z_diag_sin_primary);
%offset=0;
plot((profile_z_diag_sin_primary+1*offset),'b')
ylabel('Z (mm)')
% yt = get(gca, 'YTick');
% set(gca, 'YTick',yt, 'ZTickLabel',fliplr(yt))
%ylim([-1 1])
yyaxis right
plot(profile_t_diag)
xlabel('X (mm)')
ylabel('t (°C)')
title('Profile in x')



%quitamos perfil secundario
p1 = polyfit(1:length(profile_z_diag_sin_primary),profile_z_diag_sin_primary,3); 
f1 = polyval(p1,1:length(profile_z_diag_sin_primary)); 
figure(8);
subplot(211)
yyaxis left
plot(x_plt_dig,-(profile_z_diag_sin_primary-f1))
ylabel('Z (mm)')
%ylim([-1 1])
yyaxis right
plot(x_plt_dig,profile_t_diag)
xlabel('Profile (mm)')
ylabel('t (°C)')
%title('Digital')
grid on



 

%%
n_profile=1000;
x_linea=n_profile;
figure;

subplot(121)
plot(Y(:,n_profile),z_est(:,n_profile))
xlabel('Y (mm)')
ylabel('Z (mm)')
subplot(122)
plot(Y(:,n_profile),t_est(:,n_profile))
xlabel('Y (mm)')
ylabel('t (°C)')

% Vista de perfiles en la misma grafica sin perfil primario
profile=z_est(:,n_profile);
profile=profile(~isnan(profile));
Y_sin_nan=Y(:,n_profile);
Y_sin_nan=Y_sin_nan(~isnan(profile));


p = polyfit((1:length(profile))',profile,1); 
f = polyval(p,(1:length(profile))'); 
profile_z_sin_primary=profile-f;
% figure;
% plot(f)
% hold on
% plot(profile)

%quitamos perfil primario
%Buscamos el perfil termico
profile_t=t_est(:,n_profile);
profile_t=profile_t(~isnan(profile));

size(profile_z_sin_primary)
size(profile_t)
figure(7);
subplot(212)
yyaxis left
offset=-min(profile_z_sin_primary);
plot(Y_sin_nan,-(profile_z_sin_primary+offset))
ylabel('Z (mm)')
%ylim([-1 1])
yyaxis right
plot(Y_sin_nan,profile_t)
ylabel('t (°C)')
xlabel('Y (mm)')
title('Profile in y')


%quitamos perfil secundario
p1 = polyfit((1:length(profile_z_sin_primary))',profile_z_sin_primary,3); 
f1 = polyval(p1,(1:length(profile_z_sin_primary))'); 
figure(8);
subplot(212)
yyaxis left
plot(Y_sin_nan,-(profile_z_sin_primary-f1))
ylabel('Z (mm)')
%ylim([-1 1])
yyaxis right
plot(Y_sin_nan,profile_t)
xlabel('Y (mm)')
ylabel('t (°C)')
%%
figure;
imagesc(t_est)
colorbar
axis equal

xline(x_linea,"--k",'Profile 1','LineWidth',1)
yline(y_linea,"--k",'Profile 2','LineWidth',1)
ylim([0 1600])
c = colorbar; 
c.Label.String = '° C';
%%
figure;
yyaxis left
plot(z_sin_perfil(1000,100:end))
yyaxis right
plot(t_sin_perfil(1000,100:end))

figure;
imagesc(t_sin_perfil)
hold on
yline(1000)

%%
% % Filter sobel
% figure;
% imshow(isnan(t_est))
% z_est1=z_est(300:end-500,300:end-300);
% t_est1=t_est(300:end-500,300:end-300);
% figure; imshow(isnan(z_est))
% z_estq=histeq(z_est1);
% t_estq=histeq(t_est1);
% figure;
% subplot(121)
% imagesc(z_estq)
% subplot(121)
% imagesc(z_estq)
% 
% 
% 
% [z_est_dt,threshOutz,Gxz,Gyz] = edge(z_est1,'sobel');
% [t_est_dt,threshOutt,Gxt,Gyt]  = edge(t_est1,'sobel');
% 
% figure;
% %subplot(121)
% imagesc(Gxz)
% colorbar
% axis equal
% figure;
% imagesc(Gxt)
% colorbar
% axis equal
% 
% 
% % figure;
% % subplot(121)
% % plot(z_est_dt(1500,:))
% % subplot(122)
% % plot(t_est_dt(1500,:))
% 
% 
% 
% figure
% yyaxis left
% plot(z_est_dt(1500,:))
% yyaxis right
% plot(t_est_dt(1500,:))
% 
% 
% 
%     figure;
%            s = surf(X,Y,z_est,'FaceColor', 'interp',...
%                         'EdgeColor', 'none',...
%                         'FaceLighting', 'phong');
%     
%         
%             set(gca, 'DataAspectRatio', [1, 1, 1])
%             axis equal;
%             view(180,90);
%             camlight right
%             axis on
%             grid on
%             xlabel('x (mm)')
%             ylabel('y (mm)')
%             zlabel('z (mm)')



%% Circulos
% Digital:
aux_t_est_dig=t_est_dig;
mean_t_est=mean(mean(t_est_dig,"omitnan"),"omitnan");
aux_t_est_dig(isnan(t_est_dig))=mean_t_est;

aux_z_est_dig=z_est_dig;
mean_z_est=mean(mean(z_est_dig,"omitnan"),"omitnan");
aux_z_est_dig(isnan(z_est_dig))=mean_z_est;

sum(sum(isnan(aux_z_est_dig))) % Ojo con esto
sum(sum(isnan(aux_t_est_dig))) % Ojo con esto
save('datos_mapeo1.mat',"aux_t_est_dig","aux_z_est_dig","X_dig","Y_dig")

% Traditional
aux_t_est_trad=t_est_trad;
mean_t_est=mean(mean(t_est_trad,"omitnan"),"omitnan");
aux_t_est_trad(isnan(t_est_trad))=mean_t_est;

aux_z_est_trad=z_est_trad;
mean_z_est=mean(mean(z_est_trad,"omitnan"),"omitnan");
aux_z_est_trad(isnan(z_est_trad))=mean_z_est;

sum(sum(isnan(aux_z_est_trad))) % Ojo con esto
sum(sum(isnan(aux_t_est_trad))) % Ojo con esto
save('datos_mapeo2.mat',"aux_t_est_trad","aux_z_est_trad","X_trad","Y_trad")


%% Cargando
% Digital

%Detectando circulos en z

load([direc 'datos_mapeo_sin_forma_z_dig.mat'])
    
figure;
s = surf(XcM,YcM,Spz,'FaceColor', 'interp',...
            'EdgeColor', 'none',...
            'FaceLighting', 'phong');


set(gca, 'DataAspectRatio', [1, 1, 1])
axis equal;
view(150,45)
zt = get(gca, 'ZTick');
set(gca, 'ZTick',zt, 'ZTickLabel',fliplr(zt))
camlight left
axis on
grid on
xlabel('x (mm)')
ylabel('y (mm)')
zlabel('z (mm)')

z_sin_perfil=Spz;
figure;
imshow(z_sin_perfil,[])

sum(sum(isnan(z_sin_perfil))) % Ojo con esto

z_sin_perfil1= imadjust((z_sin_perfil));
figure;
imshow(z_sin_perfil1,[])



a=imbinarize(z_sin_perfil);
figure;
imshow(a)

figure; imagesc(a)


% figure;
% %imagesc(z_sin_perfil1)
% imshow(z_sin_perfil1,[])
% 
%  patternSize = [5, 7];  %[12, 18];         % [7, 21]
%  Pattern_type = 'asymmetric';  %'symmetric';     % asymmetric
%  Circ_color = 'black';           % white
%  puntos = detectCircleGridPoints(uint8(a),patternSize,PatternType=Pattern_type,CircleColor=Circ_color);
% figure;
% imagesc(z_sin_perfil)
% hold on
% plot(puntos(:,1),puntos(:,2),'+r')

[centers, radii, metric] = imfindcircles(z_sin_perfil,[40 120],'ObjectPolarity','dark',Sensitivity=0.928);
figure;
imagesc(a)
hold on
plot(centers(:,1),centers(:,2),'ok')
viscircles(centers,radii)


% aux_z_est=z_est;
% mean_z_est=mean(mean(z_est,"omitnan"),"omitnan");
% aux_z_est(isnan(z_est))=mean_z_est;
% %aux_z_est=histeq(aux_z_est);
% figure;
% imagesc(aux_z_est)

figure; imagesc(z_sin_perfil1)



%%
% figure;
% subplot(121)
% imagesc(a)
% axis equal
% hold on
% subplot(122)
% imagesc(aux_t_est)
% axis equal
% hold on
% for i=1:size(centers,1)
%     subplot(121)
%     plot(centers(i,1),centers(i,2),'ok')
%     subplot(122)
%     plot(puntos(i,1),puntos(i,2),'+r')
%     pause();
% end



%%
centers1=centers;
% Buscamos el inicio
centers_ord=[];
% for i=1:size(centers1,1)
%     distance(i)=norm(centers1(i,:)- [0 size(a,1)]);
% end
distance=vecnorm((centers1-[0 size(a,1)])');
[first, idx]=min(distance);
centers_ord(1,:)=centers1(idx,:);
% Buscamos los colineales horizontalmente
distance_row=abs(centers1(:,2)-centers_ord(1,2))
[B,idxs] = sort(distance_row);
candidates=centers1(idxs(1:5),:)
centers1(idxs(1:5),:)=[]; % Borramos los candidates del listado de centros
search_idx=B==0;
candidates(search_idx(1:5),:)=[];
% Ordenamos los candidatos
distance_candidates=[];
% for i=1:size(candidates,1)
%     distance_candidates(i)=norm(candidates(i,:)- centers_ord(1,:));
% end
distance_candidates=vecnorm((candidates-centers_ord(1,:))');
[B,idxs] = sort(distance_candidates);
centers_ord(size(centers_ord,1)+1:size(centers_ord,1)+4,:)=candidates(idxs,:);


for j=2:7
% Buscamos la otra fila cercana
    distance=[];
    distance=vecnorm((centers1-centers_ord(1,:))');
%     for i=1:size(centers1,1)
%         distance(i)=norm(centers1(i,:)- centers_ord(1,:));
%     end
    [first, idx]=min(distance);
    second=centers1(idx,:);
    % figure;
    % imagesc(aux_t_est)
    % hold on
    % plot(second(1,1), second(1,2),'ok')
    
    centers_ord(size(centers_ord,1)+1,:)=centers1(idx,:);
    % Buscamos los colineales horizontalmente
    distance_row=abs(centers1(:,2)-centers_ord(size(centers_ord,1),2))
    [B,idxs] = sort(distance_row);
    candidates=centers1(idxs(1:5),:)
    centers1(idxs(1:5),:)=[]; % Borramos los candidates del listado de centros
    search_idx=B==0;
    candidates(search_idx(1:5),:)=[];
    % Ordenamos los candidatos
    distance_candidates=[];
    distance_candidates=vecnorm((candidates-centers_ord(1,:))');
%     for i=1:size(candidates,1)
%         distance_candidates(i)=norm(candidates(i,:)- centers_ord(1,:));
%     end
    [B,idxs] = sort(distance_candidates);
    centers_ord(size(centers_ord,1)+1:size(centers_ord,1)+4,:)=candidates(idxs,:);
    

end
figure;
imagesc(aux_t_est)
hold on
for i=1:size(centers_ord,1)
    plot(centers_ord(i,1),centers_ord(i,2),'ok')
    pause(0.2);
end
puntos_z=centers_ord;
% %%
% for j=3:7
%     % Buscamos la tercera fila
%     distance=[];
%     for i=1:size(centers1,1)
%         distance(i)=norm(centers1(i,:)- centers_ord(1,:));
%     end
%     [first, idx]=min(distance);
%     third=centers1(idx,:);
%     figure;
%     imagesc(aux_t_est)
%     hold on
%     plot(third(1,1), third(1,2),'ok')
%     centers_ord(size(centers_ord,1)+1,:)=centers1(idx,:);
%     % Buscamos los colineales horizontalmente
%     distance_row=abs(centers1(:,2)-centers_ord(size(centers_ord,1),2))
%     [B,idxs] = sort(distance_row);
%     candidates=centers1(idxs(1:5),:)
%     centers1(idxs(1:5),:)=[]; % Borramos los candidates del listado de centros
%     search_idx=B==0;
%     candidates(search_idx(1:5),:)=[];
%     % Ordenamos los candidatos
%     distance_candidates=[];
%     for i=1:size(candidates,1)
%         distance_candidates(i)=norm(candidates(i,:)- centers_ord(1,:));
%     end
%     [B,idxs] = sort(distance_candidates);
%     centers_ord(size(centers_ord,1)+1:size(centers_ord,1)+4,:)=candidates(idxs,:);
%     
%     figure;
%     imagesc(aux_t_est)
%     hold on
%     for i=1:size(centers_ord,1)
%         plot(centers_ord(i,1),centers_ord(i,2),'ok')
%         pause();
%     end
% end

%%
% centers_ord = [];
% % Buscamos el inicio
% distances = vecnorm(centers - [0 size(a,1)], 2, 2);
% [~, idx] = min(distances);
% centers_ord(1,:) = centers(idx,:);
% centers(idx,:) = []; % Eliminamos el centro encontrado
% 
% % Función auxiliar para ordenar candidatos
% ordenar_candidatos = @(cands, ref) sortrows(cands, vecnorm(cands - ref, 2, 2));
% 
% % Procesamiento de filas
% for j = 1:7
%     % Buscamos la fila cercana
%     distances = vecnorm(centers - centers_ord(1,:), 2, 2);
%     [~, idx] = min(distances);
%     fila_cercana = centers(idx,:);
%     centers_ord(end+1,:) = fila_cercana;
%     centers(idx,:) = []; % Eliminamos el centro encontrado
% 
%     % Buscamos los colineales horizontalmente
%     tol=0.00001
%     idxs_fila = abs(centers(:,2) - fila_cercana(2)) < tol; % Ajusta tol según sea necesario
%     candidatos = centers(idxs_fila,:);
%     centers(idxs_fila,:) = []; % Eliminamos los candidatos del listado de centros
% 
%     % Ordenamos los candidatos y agregamos
%     candidatos_ordenados = ordenar_candidatos(candidatos, fila_cercana);
%     centers_ord = [centers_ord; candidatos_ordenados(1:4,:)]; % Añadir los 4 más cercanos
%     
%     % Visualización
%     figure;
%     imagesc(aux_t_est)
%     hold on
%     plot(centers_ord(:,1), centers_ord(:,2), 'ok')
%     pause();
% end


%%

load(["datos_mapeo_sin_forma_t.mat"])
    
figure;
s = surf(XcM,YcM,Spz,'FaceColor', 'interp',...
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

t_sin_perfil=Spz;
figure; 
imshow(t_sin_perfil,[])

t_sin_perfil1= imadjust(t_sin_perfil);
a=(imbinarize(t_sin_perfil1));
figure; 
imshow(a,[])
[puntos, radii, metric] = imfindcircles(a,[84 140],'ObjectPolarity','bright',Sensitivity=0.96);
figure;
imagesc(a)
hold on
plot(puntos(:,1),puntos(:,2),'ok')
viscircles(puntos,radii)

% patternSize = [5, 7];  %[12, 18];         % [7, 21]
% centerDistance = 40;            % 20
% Pattern_type = 'asymmetric';  %'symmetric';     % asymmetric
% Circ_color = 'black';           % white
% puntos = detectCircleGridPoints(uint8(a),patternSize,PatternType=Pattern_type,CircleColor=Circ_color);
% figure;
% imagesc(t_sin_perfil1)
% hold on
% plot(puntos(:,1),puntos(:,2),'+r')

%%
aux_z_est=z_est(20:end-20,20:end-20);
x1=X(20:end-20,20:end-20);
y1=Y(20:end-20,20:end-20);
%mean_z_est=mean(mean(z_est,"omitnan"),"omitnan");
%aux_z_est(isnan(z_est))=mean_z_est;
tam=size(aux_z_est);
x2=x1(:);
y2=y1(:);
z1=aux_z_est(:);
DM = [x2, y2, ones(size(z1))];                             % Design Matrix
%B = DM\z;                                               % Estimate Parameters
B=lsqminnorm(DM,z1);
%[X,Y] = meshgrid(linspace(min(x1),max(x1),50), linspace(min(y),max(y),50));
Z = B(1)*x1 + B(2)*y1 + B(3)*ones(size(x1));

z_sin_forma=aux_z_est-Z;

z_sin_forma(z_sin_forma>max(z))=NaN;
figure;
s = surf(x1,y1,z_sin_forma,'FaceColor', 'interp',...
            'EdgeColor', 'none',...
            'FaceLighting', 'phong');


set(gca, 'DataAspectRatio', [1, 1, 1])
axis equal;
view(150,45)
zt = get(gca, 'ZTick');
set(gca, 'ZTick',zt, 'ZTickLabel',fliplr(zt))
camlight left
axis on
grid on
xlabel('x (mm)')
ylabel('y (mm)')
zlabel('z (mm)')

figure;
imagesc(z_sin_forma)
hold on
yline(970,'--k')

c = colorbar; 
c.Label.String = 'z (mm)';

figure;
plot(x1(970,:),-z_sin_forma(970,:))
xlabel('x (mm)')
ylabel('z (mm)')
xlim([-110 100])
grid on
%yt = get(gca, 'YTick');
%set(gca, 'YTick',yt, 'YTickLabel',fliplr(yt))

%%

centers1=puntos;
% Buscamos el inicio
centers_ord=[];
% for i=1:size(centers1,1)
%     distance(i)=norm(centers1(i,:)- [0 size(a,1)]);
% end
distance=vecnorm((centers1-[0 size(a,1)])');
[first, idx]=min(distance);
centers_ord(1,:)=centers1(idx,:);
% Buscamos los colineales horizontalmente
distance_row=abs(centers1(:,2)-centers_ord(1,2))
[B,idxs] = sort(distance_row);
candidates=centers1(idxs(1:5),:)
centers1(idxs(1:5),:)=[]; % Borramos los candidates del listado de centros
search_idx=B==0;
candidates(search_idx(1:5),:)=[];
% Ordenamos los candidatos
distance_candidates=[];
% for i=1:size(candidates,1)
%     distance_candidates(i)=norm(candidates(i,:)- centers_ord(1,:));
% end
distance_candidates=vecnorm((candidates-centers_ord(1,:))');
[B,idxs] = sort(distance_candidates);
centers_ord(size(centers_ord,1)+1:size(centers_ord,1)+4,:)=candidates(idxs,:);


for j=2:7
% Buscamos la otra fila cercana
    distance=[];
    distance=vecnorm((centers1-centers_ord(1,:))');
%     for i=1:size(centers1,1)
%         distance(i)=norm(centers1(i,:)- centers_ord(1,:));
%     end
    [first, idx]=min(distance);
    second=centers1(idx,:);
    % figure;
    % imagesc(aux_t_est)
    % hold on
    % plot(second(1,1), second(1,2),'ok')
    
    centers_ord(size(centers_ord,1)+1,:)=centers1(idx,:);
    % Buscamos los colineales horizontalmente
    distance_row=abs(centers1(:,2)-centers_ord(size(centers_ord,1),2))
    [B,idxs] = sort(distance_row);
    candidates=centers1(idxs(1:5),:)
    centers1(idxs(1:5),:)=[]; % Borramos los candidates del listado de centros
    search_idx=B==0;
    candidates(search_idx(1:5),:)=[];
    % Ordenamos los candidatos
    distance_candidates=[];
    distance_candidates=vecnorm((candidates-centers_ord(1,:))');
%     for i=1:size(candidates,1)
%         distance_candidates(i)=norm(candidates(i,:)- centers_ord(1,:));
%     end
    [B,idxs] = sort(distance_candidates);
    centers_ord(size(centers_ord,1)+1:size(centers_ord,1)+4,:)=candidates(idxs,:);
    

end
figure;
imagesc(t_sin_perfil)
hold on
for i=1:size(centers_ord,1)
    plot(centers_ord(i,1),centers_ord(i,2),'ok')
    pause(0.2);
end
puntos_t=centers_ord;


%%

z_est_trad_without_surface=[];
t_est_trad_without_surface=[];
for i=1:size(z_est_trad,1)
    % Digital Features
    %close all
    n_profile=i;
    
    largo_profile=length(t_est_trad(i,:));
    figure(1);
    imagesc(t_est_trad)
    yline(n_profile,'--k')
    colorbar
    
    
    tam=length(X_trad(n_profile,:));

    y_linea=n_profile;
%     figure;
%     subplot(121)
%     plot(z_est_trad(n_profile,largo_profile))
%     xlabel('X (mm)')
%     ylabel('Z (mm)')
%     subplot(122)
%     plot(t_est_trad(n_profile,largo_profile))
%     xlabel('X (mm)')
%     ylabel('t (°C)')
    
    
    % Vista de perfiles en la misma grafica sin perfil primario
    profile_real=-z_est_trad(n_profile,120:tam-20);
    mask=isnan(profile_real);
    profile=profile_real(~mask);
    X_trad_sin_nan=X_trad(n_profile,:);
    X_trad_sin_nan=X_trad_sin_nan(~mask);
    
    
    p = polyfit(1:length(profile),profile,1); 
    f = polyval(p,1:length(profile)); 
    profile_z_trad_sin_primary=profile-f;
    % figure;
    % plot(f)
    % hold on
    % plot(profile)
    
    %quitamos perfil primario
    %Buscamos el perfil termico
    profile_t_trad=t_est_trad(n_profile,120:tam-20);
    profile_t_trad=profile_t_trad(~mask);
    
    size(profile_z_trad_sin_primary)
    size(profile_t_trad)
    figure(7);
    subplot(212)
    yyaxis left
    offset=-min(profile_z_trad_sin_primary);
    %offset=0;
    plot((profile_z_trad_sin_primary+1*offset),'b')
    ylabel('Z (mm)')
    % yt = get(gca, 'YTick');
    % set(gca, 'YTick',yt, 'ZTickLabel',fliplr(yt))
    %ylim([-1 1])
    yyaxis right
    plot(profile_t_trad)
    xlabel('X (mm)')
    ylabel('t (°C)')
    title('Profile in x')
    
    
    
    %quitamos perfil secundario
    p1 = polyfit(1:length(profile_z_trad_sin_primary),profile_z_trad_sin_primary,3); 
    f2 = polyval(p1,1:length(profile_z_trad_sin_primary)); 
    figure(8);
    subplot(212)
    yyaxis left
    plot((profile_z_trad_sin_primary-f2))
    ylabel('Z (mm)')
    %ylim([-1 1])
    yyaxis right
    plot(profile_t_trad)
    xlabel('X (mm)')
    ylabel('t (°C)')
    title('Traditional')
    grid on
    if length(profile_z_trad_sin_primary-f2)>0
        z_est_trad_without_surface(i,~mask)=(profile_z_trad_sin_primary-f2);
        t_est_trad_without_surface(i,~mask)=profile_t_trad;
    else
        z_est_trad_without_surface(i,:)=NaN*ones(1,length(largo_profile));
        t_est_trad_without_surface(i,:)=NaN*ones(1,length(largo_profile));
    end
    pause(0.1)
    %clf(7)
    %clf(8)
end

figure; imagesc(z_est_trad_without_surface)
figure; imagesc(t_est_trad_without_surface)


%%
close all
t_sin_perfil=t_est_trad_without_surface;
figure; 
imshow(t_sin_perfil,[])
B = sort(t_sin_perfil(:));
C = unique(B);
t_sin_perfil(isnan(t_sin_perfil))=C(2); %mean(mean(t_sin_perfil,"omitnan"),"omitnan");
t_sin_perfil(t_sin_perfil==0)=C(2); %mean(mean(t_sin_perfil,"omitnan"),"omitnan");

%t_sin_perfil1= imadjust(t_sin_perfil,[0.99 1]);
figure; 
imshow(t_sin_perfil,[])
%t_sin_perfil1= (t_sin_perfil);

X=t_sin_perfil;
Xmin = min(X(:));
Xmax = max(X(:));
if isequal(Xmax,Xmin)
    X = 0*X;
else
    X = (X - Xmin) ./ (Xmax - Xmin);
end

% Threshold image - adaptive threshold
BW = imbinarize(X, 'adaptive', 'Sensitivity', 0.500000, 'ForegroundPolarity', 'bright');

% Open mask with disk
radius = 10;
decomposition = 0;
se = strel('disk', radius, decomposition);
a = imopen(BW, se);

% Threshold image - adaptive threshold
figure; imshow(a)
[mask_w_w,numWhite] = bwlabel(a);
mask_w_w = bwareaopen(mask_w_w,10000);  %Retirando manchas blancas
figure; imshow(mask_w_w)

%a=(imbinarize(t_sin_perfil));
a=mask_w_w ;
% figure; 
% imshow(a,[])
% [puntos, radii, metric] = imfindcircles(a,[50 70],'ObjectPolarity','bright',Sensitivity=0.9);
% figure;
% imagesc(a)
% hold on
% plot(puntos(:,1),puntos(:,2),'ok')
% viscircles(puntos,radii)

figure;
imagesc(a)
stats = regionprops("table",a,"Centroid", ...
    "MajorAxisLength","MinorAxisLength")
centers_t = stats.Centroid;
diameters = mean([stats.MajorAxisLength stats.MinorAxisLength],2);
radii = diameters/2;
hold on
viscircles(centers_t,radii)
plot(centers_t(:,1),centers_t(:,2),'+k')


%%
close all
z_sin_perfil=z_est_trad_without_surface;
figure; 
imshow(z_sin_perfil,[])
B = sort(z_sin_perfil(:));
C = unique(B);
z_sin_perfil(isnan(z_sin_perfil))=C(2); %mean(mean(t_sin_perfil,"omitnan"),"omitnan");
z_sin_perfil(z_sin_perfil==0)=C(2); %mean(mean(t_sin_perfil,"omitnan"),"omitnan");


X=imsharpen(z_sin_perfil);
%X=imadjust(z_sin_perfil);
Xmin = min(X(:));
Xmax = max(X(:));
if isequal(Xmax,Xmin)
    X = 0*X;
else
    X = (X - Xmin) ./ (Xmax - Xmin);
end

% BW1 = edge(X,"approxcanny");
% figure; 
% imshow(BW1,[])
% se = strel('disk',3);
% closeBW = imclose(BW1,se);
% figure, imshow(closeBW)
% 
% BW2 = imfill(closeBW,8,"holes"); figure,    imshow(BW2)
% 
%     imshow(BW2)



X=imadjust(z_sin_perfil);
% Threshold image - adaptive threshold
BW = imbinarize(X, 'adaptive', 'Sensitivity', 0.700000, 'ForegroundPolarity', 'bright');
figure; 
imshow(BW,[])
BW = imfill(BW,"holes"); figure,    imshow(BW)

[mask_w_w,numWhite] = bwlabel(BW);
mask_w_w = bwareaopen(mask_w_w,7000);  %Retirando manchas blancas
figure; imshow(mask_w_w)
[mask_w_w,numWhite] = bwlabel(mask_w_w);

stats = regionprops("table",mask_w_w,"Centroid", ...
    "MajorAxisLength","MinorAxisLength",'Area')
centers_t = stats.Area;
[r, c] = find(L==2);
rc = [r c]

[puntos, radii, metric] = imfindcircles(BW,[50 70],'ObjectPolarity','bright',Sensitivity=0.9);
figure;
imagesc(BW)
hold on
plot(puntos(:,1),puntos(:,2),'ok')
viscircles(puntos,radii)
%%
z_sin_perfil=imsharpen(z_sin_perfil);
z_sin_perfil=imadjust(z_sin_perfil);
%t_sin_perfil1= imadjust(t_sin_perfil,[0.99 1]);
figure; 
imshow(z_sin_perfil,[])
%t_sin_perfil1= (t_sin_perfil);

X=z_sin_perfil;
Xmin = min(X(:));
Xmax = max(X(:));
if isequal(Xmax,Xmin)
    X = 0*X;
else
    X = (X - Xmin) ./ (Xmax - Xmin);
end

% Threshold image - adaptive threshold
BW = imbinarize(X, 'adaptive', 'Sensitivity', 0.500000, 'ForegroundPolarity', 'bright');

% Open mask with disk
radius = 10;
decomposition = 0;
se = strel('disk', radius, decomposition);
a = imopen(BW, se);

% Threshold image - adaptive threshold
figure; imshow(a)
[mask_w_w,numWhite] = bwlabel(a);
mask_w_w = bwareaopen(mask_w_w,10000);  %Retirando manchas blancas
figure; imshow(mask_w_w)

%a=(imbinarize(t_sin_perfil));
a=mask_w_w ;
% figure; 
% imshow(a,[])
% [puntos, radii, metric] = imfindcircles(a,[50 70],'ObjectPolarity','bright',Sensitivity=0.9);
% figure;
% imagesc(a)
% hold on
% plot(puntos(:,1),puntos(:,2),'ok')
% viscircles(puntos,radii)

figure;
imagesc(a)
stats = regionprops("table",a,"Centroid", ...
    "MajorAxisLength","MinorAxisLength")
centers_t = stats.Centroid;
diameters = mean([stats.MajorAxisLength stats.MinorAxisLength],2);
radii = diameters/2;
hold on
viscircles(centers_t,radii)
plot(centers_t(:,1),centers_t(:,2),'+k')
