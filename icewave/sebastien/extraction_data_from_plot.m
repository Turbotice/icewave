%% This script decprits how we can extract data from a matlab figure

ax = gca; 
h = findobj(gca,'Type','line');

fitted_alpha = h(2).YData;
fitted_f = h(2).XData;

f_list = h(1).XData;
yth = h(1).YData;

%% Save variables 
filename = 'Data_attenuation_law_distance_0p05';
fig_folder = 'W:/Banquise/Rimouski_2023/Traitements_donnees/PIV_Sebastien/Figures_report/';

save([fig_folder filename],'fitted_alpha','fitted_f','f_list','yth','-v7.3')

%% Check content of some variables 

ax = gca; 
h = findobj(gca,'Type','line');

%%
figure,

plot(h(2).XData,h(2).YData)
hold on 
plot(h(1).XData,h(1).YData,'.')
xlabel('$y \: \rm (cm)$')
ylabel('$h \: \rm (cm)$')
ax = gca; 
ax.FontSize = 13;
set_Papermode(gcf);

