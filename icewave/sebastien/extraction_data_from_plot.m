%% This script decprits how we can extract data from a matlab figure


%% Extract data from Figure 5f of article 
ax = gca; 
h = findobj(gca,'Type','line');

envelop = struct();
envelop.x = h(1).XData;
envelop.y = h(1).YData;

demod = struct();
demod.x = h(2).XData;
demod.y = h(2).YData;

profile_y = struct();
profile_y.envelop = envelop;
profile_y.demod = demod;

base = 'C:/Users/sebas/OneDrive/Bureau/These PMMH/Waves_float/Figures_article/';
fig_folder =  [base 'Inkscape_drawings/figure5_referees_corrected/'];
filename = [fig_folder 'profile_Y_data'];

save(filename,'profile_y','-v7.3')

%% Extract data from X profile - figure 5e of article
ax = gca; 
h = findobj(gca,'Type','line');

envelop = struct();
envelop.x = h(1).XData;
envelop.y = h(1).YData;

demod = struct();
demod.x = h(2).XData;
demod.y = h(2).YData;

profile_x = struct();
profile_x.envelop = envelop;
profile_x.demod = demod;

base = 'C:/Users/sebas/OneDrive/Bureau/These PMMH/Waves_float/Figures_article/';
fig_folder =  [base 'Inkscape_drawings/figure5_referees_corrected/'];
filename = [fig_folder 'profile_X_data'];

save(filename,'profile_x','-v7.3')




%%
figure, 
plot(demod.x,demod.y)
hold on 
plot(envelop.x,envelop.y,'.')






%% Previous extraction
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

