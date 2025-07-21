function [z_est_final,t_est_final]=detrend_profile_y(z_est, t_est)





tam=length(t_est);
%tam=length(X_trad(n_profile,:));

% Vista de perfiles en la misma grafica sin perfil primario
profile_real=-z_est;
mask=isnan(profile_real);
profile=profile_real(~mask);



p = polyfit(1:length(profile),profile',2); 
f = polyval(p,1:length(profile)); 
profile_z_sin_primary=profile-f';
% figure;
% plot(f)
% hold on
% plot(profile)

%quitamos perfil primario
%Buscamos el perfil termico
profile_t=t_est;
profile_t=profile_t(~mask);

size(profile_z_sin_primary);
size(profile_t);
% figure(7);
% subplot(212)
% yyaxis left
% offset=-min(profile_z_sin_primary);
% offset=0;
% plot((profile_z_sin_primary+1*offset),'b')
% ylabel('Z (mm)')
% % yt = get(gca, 'YTick');
% % set(gca, 'YTick',yt, 'ZTickLabel',fliplr(yt))
% %ylim([-1 1])
% yyaxis right
% plot(profile_t)
% xlabel('X (mm)')
% ylabel('t (°C)')
% title('Profile in x')



%quitamos perfil secundario
p1 = polyfit(1:length(profile_z_sin_primary),profile_z_sin_primary',3); 
f2 = polyval(p1,1:length(profile_z_sin_primary)); 
profile_z_final=profile_z_sin_primary-f2';
% figure(8);
% subplot(212)
% yyaxis left
% plot((profile_z_sin_primary-f2))
% ylabel('Z (mm)')
% %ylim([-1 1])
% yyaxis right
% plot(profile_t)
% xlabel('X (mm)')
% ylabel('t (°C)')
% title('Traditional')
% grid on
% if ~isempty(profile_z_sin_primary-f2)
%     z_est_final=NaN*ones(1,tam);
%     t_est_final=NaN*ones(1,tam);
%     z_est_final(~mask)=(profile_z_sin_primary-f2);
%     t_est_final(~mask)=profile_t;
% else
%     z_est_final=NaN*ones(1,tam);
%     t_est_final=NaN*ones(1,tam);
% end

if ~isempty(profile_z_sin_primary)
    z_est_final=NaN*ones(1,tam);
    t_est_final=NaN*ones(1,tam);
    z_est_final(~mask)=(profile_z_sin_primary);
    t_est_final(~mask)=profile_t;
else
    z_est_final=NaN*ones(1,tam);
    t_est_final=NaN*ones(1,tam);
end










end