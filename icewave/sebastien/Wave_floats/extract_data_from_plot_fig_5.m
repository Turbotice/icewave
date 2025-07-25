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





