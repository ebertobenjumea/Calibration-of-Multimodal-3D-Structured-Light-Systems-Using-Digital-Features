close all
clc
clear all
direc='D:\R3D\Adquisiciones\Experimento 54\Objects\SL\';
aux=num2str(35);
% Cargamos digital features image
pc_dig = pcread([direc 'Digital Features 1 UTB\new_obj_vis_ir_' aux '.ply']);
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
pc_trad = pcread([direc 'Two Stages 2\new_obj_vis_ir_' aux '.ply']);
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


%%  digital
close all
clc
% Suponiendo que tienes los puntos como vectores Nx1
% X, Y, Z

% Combina en una sola matriz
data = [x_dig, y_dig, -z_dig];
mean_data=mean(data);
data=data-mean_data;

% Estimación inicial del centro y radio
x0 = mean(data);
r0 = mean(sqrt(sum((data - x0).^2, 2)));
params0 = [x0, r0];  % [xc, yc, zc, r]

% Función de error: distancia desde el punto al centro menos el radio
error_fun = @(params) sqrt(sum((data - params(1:3)).^2, 2)) - params(4);

% Ajuste por mínimos cuadrados no lineales
options = optimoptions('lsqnonlin', 'Display', 'off');
params_fit = lsqnonlin(error_fun, params0, [], [], options);

% Resultado
center = params_fit(1:3);
radius = params_fit(4);


% Mostrar puntos originales
figure;
scatter3(data(:,1), data(:,2), data(:,3),  ...
        'MarkerFaceColor',[0 0 1],'MarkerFaceAlpha',0.5);
hold on;

% Dibujar la esfera ajustada
[xs, ys, zs] = sphere(50);
surf(radius*xs + center(1), radius*ys + center(2), radius*zs + center(3), ...
    'FaceAlpha', 1, 'EdgeColor', 'none',FaceColor="yellow");

axis equal;
grid on;
%title('Ajuste de esfera a nube de puntos');


%%






% Aplicar interpolación de colores y sin líneas de malla
shading interp

% Añadir una fuente de luz
light

% Configurar el tipo de iluminación
lighting gouraud  % O prueba 'phong' para efectos más suaves

% Configurar el material para controlar cómo refleja la luz
material shiny     % También puedes probar 'dull', 'metal', etc.

% Otras mejoras visuales
%axis tight
%camlight left  % Posicionar la luz desde la izquierda
camlight(0,270)
%camlight("left","infinite")


% Estimación inicial del centro y radio
x0 = mean(data);
r0 = mean(sqrt(sum((data - x0).^2, 2)));
params0 = [x0, r0];  % [xc, yc, zc, r]

% Función de error: diferencia entre distancia al centro y el radio
error_fun = @(params) sqrt(sum((data - params(1:3)).^2, 2)) - params(4);

% Ajuste con lsqnonlin
options = optimoptions('lsqnonlin', 'Display', 'off');
params_fit = lsqnonlin(error_fun, params0, [], [], options);

% Resultado del ajuste
center = params_fit(1:3);
radius = params_fit(4);

% --- CÁLCULO DEL ERROR RMSE ---
dist_to_center = sqrt(sum((data - center).^2, 2));
errors = dist_to_center - radius;
rmse = sqrt(mean(errors.^2));
figure;
histogram(errors,"EdgeColor","none")
xlabel('Error (mm)')
ylabel('Counts')

% Mostrar resultado
fprintf('Centro ajustado: (%.4f, %.4f, %.4f)\n', center);
fprintf('Radio ajustado: %.4f\n', radius);
fprintf('RMSE del ajuste: %.6f\n', rmse);

%% Tradicional

% Suponiendo que tienes los puntos como vectores Nx1
% X, Y, Z

% Combina en una sola matriz
data = [x_trad, y_trad, -z_trad];
mean_data=mean(data);
data=data-mean_data;
% Estimación inicial del centro y radio
x0 = mean(data);
r0 = mean(sqrt(sum((data - x0).^2, 2)));
params0 = [x0, r0];  % [xc, yc, zc, r]

% Función de error: distancia desde el punto al centro menos el radio
error_fun = @(params) sqrt(sum((data - params(1:3)).^2, 2)) - params(4);

% Ajuste por mínimos cuadrados no lineales
options = optimoptions('lsqnonlin', 'Display', 'off');
params_fit = lsqnonlin(error_fun, params0, [], [], options);

% Resultado
center = params_fit(1:3);
radius = params_fit(4);


% Mostrar puntos originales
%figure;
scatter3(data(:,1), data(:,2), data(:,3), ...
         'MarkerFaceColor',[1 0 0],'MarkerFaceAlpha',1);
%hold on;

% Dibujar la esfera ajustada
[xs, ys, zs] = sphere(50);
%surf(radius*xs + center(1), radius*ys + center(2), radius*zs + center(3), ...
%    'FaceAlpha', 0.3, 'EdgeColor', 'none', 'FaceColor', 'r');

axis equal;
grid on;
%title('Ajuste de esfera a nube de puntos');
xlabel('x (mm)')
ylabel('y (mm)')
zlabel('z (mm)')



% Estimación inicial del centro y radio
x0 = mean(data);
r0 = mean(sqrt(sum((data - x0).^2, 2)));
params0 = [x0, r0];  % [xc, yc, zc, r]

% Función de error: diferencia entre distancia al centro y el radio
error_fun = @(params) sqrt(sum((data - params(1:3)).^2, 2)) - params(4);

% Ajuste con lsqnonlin
options = optimoptions('lsqnonlin', 'Display', 'off');
params_fit = lsqnonlin(error_fun, params0, [], [], options);

% Resultado del ajuste
center = params_fit(1:3);
radius = params_fit(4);

% --- CÁLCULO DEL ERROR RMSE ---
dist_to_center = sqrt(sum((data - center).^2, 2));
errors = dist_to_center - radius;
rmse = sqrt(mean(errors.^2));

% Mostrar resultado
fprintf('Centro ajustado: (%.4f, %.4f, %.4f)\n', center);
fprintf('Radio ajustado: %.4f\n', radius);
fprintf('RMSE del ajuste: %.6f\n', rmse);

%camlight(90,0)
legend('Digital features','Fitted sphere','Conventional')