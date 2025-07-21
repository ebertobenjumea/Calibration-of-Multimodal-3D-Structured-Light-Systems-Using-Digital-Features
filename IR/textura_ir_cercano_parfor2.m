function [x_hat, y_hat,x,y,z_hat, ir_texture, nan_idx]=textura_ir_cercano_parfor2(x_hat,y_hat,x,y,z_hat,ir_texture,img_ir)
% x_hat: points in x image coordinate of the thermal camera. (1xn)
% y_hat: points in y image coordinate of the thermal camera. (1xn)
% z_hat: z coordinate of the 3D data in thermal camera system. (1xn)
% ir_texture: temperature (1xn)
%a=1
nan_idx=zeros(size(x_hat));
nan_idx=nan_idx | isnan(x_hat);
nan_idx=nan_idx | y_hat>192;
nan_idx=nan_idx | x_hat>256;
nan_idx=nan_idx | y_hat<1;
nan_idx=nan_idx | x_hat<1;


x_hat_copy=x_hat;
x_hat=x_hat(~isnan(x_hat_copy));
y_hat=y_hat(~isnan(x_hat_copy));
x=x(~isnan(x_hat_copy));
y=y(~isnan(x_hat_copy));
z_hat=z_hat(~isnan(x_hat_copy));
ir_texture=ir_texture(~isnan(x_hat_copy));




y_hat_copy=y_hat;

x(y_hat>192)=[];
y(y_hat>192)=[];
z_hat(y_hat>192)=[];
ir_texture(y_hat>192)=[];
x_hat(y_hat>192)=[];
y_hat(y_hat>192)=[];


y_hat(x_hat>256)=[];
x(x_hat>256)=[];
y(x_hat>256)=[];
z_hat(x_hat>256)=[];
ir_texture(x_hat>256)=[];
x_hat(x_hat>256)=[];


x(y_hat<1)=[];
y(y_hat<1)=[];
z_hat(y_hat<1)=[];
ir_texture(y_hat<1)=[];
x_hat(y_hat<1)=[];
y_hat(y_hat<1)=[];

y_hat(x_hat<1)=[];
x(x_hat<1)=[];
y(x_hat<1)=[];
z_hat(x_hat<1)=[];
ir_texture(x_hat<1)=[];
x_hat(x_hat<1)=[];
%
temp = img_ir; % IR image
[m, n] = size(temp);
%m=1600;
%n=2048;
%y_hat = [1   2.3 2.1 1.5 3.1 4.5 2.1 1.6];
%x_hat = [1.1 2.2 2.5 2   2.8 1.4 2.9 1.6];
%z_hat = [100 121 120 120 119 123 130 137];
%ir_texture =[1 2 2 4 3 4 5 6];

x_r = floor(x_hat);
y_r = floor(y_hat);

ind = sub2ind([m n],y_r,x_r);





idx_erase=zeros(1,length(z_hat));
idx_erase_correct1=idx_erase;
%tic
parfor k = 1:numel(temp)
        
    if  mod(k,1000)==0
        disp([num2str(k*100/numel(temp)) '%'])
    end
    if sum(k == ind)>1
        %k
        k == ind;
        ind_z = 1:numel(z_hat);
        aux = ind_z(k == ind);
        [zmin_val, ind_zmin] = min(z_hat(k == ind));
        idx_erase=idx_erase+(k == ind);
        idx=find(~(z_hat(k == ind) > zmin_val)==1);
        %idx_erase(aux(idx))=0;
        idx_erase_correct=[zeros(1,aux(idx)-1) -1 zeros(1,length(z_hat)-aux(idx))];
        idx_erase_correct1=idx_erase_correct+idx_erase_correct1;
        %idx_erase1=idx_erase+idx_erase_correct;
        %idx_erase
        % aux(z_hat(k == ind) > zmin_val)
        %z_hat
    end
end
disp('Por aca')
idx_erase1=idx_erase+idx_erase_correct1;
idx_erase1=~logical(idx_erase1);
cnt_valid=0;
% for m=1:length(nan_idx)
%     if nan_idx(m)==0
%         cnt_valid=cnt_valid+1;
%         nan_idx(m)=~idx_erase1(cnt_valid);
%     end
% end
% nan_idx  : lógico, tamaño = nube original (true = punto inválido)
% idx_erase1 : lógico, tamaño = puntos válidos (true = conservar)

% --- Comprobación opcional de coherencia -------------------------------
%assert(nnz(~nan_idx) == numel(idx_erase1), ...
%       'idx_erase1 debe tener la misma longitud que los puntos válidos');

% --- Sustituir en un solo paso -----------------------------------------
nan_idx(~nan_idx) = ~idx_erase1;   % true  → descartar
                                   % false → conservar



x_hat=x_hat(idx_erase1);
y_hat=y_hat(idx_erase1);
z_hat=z_hat(idx_erase1);
ir_texture=ir_texture(idx_erase1);
x=x(idx_erase1);
y=y(idx_erase1);
%toc
