%% Load data 

filename = 'Data_attenuation_merged_0p1Hz_to_0p8Hz_mode_1.mat';
base = 'W:/SagWin2024/Data/0226/Drones/mesange/matData/10-waves_005/Plots/';

S = load([base filename]);

filename = 'Data_attenuation_law_distance_0p05.mat';
base = 'W:/Banquise/Rimouski_2023/Traitements_donnees/PIV_Sebastien/Figures_report/DJI_0402/';

S2 = load([base filename]);

%% mask datas 
%  masking according to dist_fit
d_thresh = 0.05;

mask = (s.dist < d_thresh) & (abs(s.alpha) > 0.001) ; % keep only points for which distance is smaller than..
f = s.freq_array;
S.fitted_f = f(mask);
S.fitted_alpha = abs(s.alpha(mask))';

%%
l1(1,:) = fminsearch(@(s)powerfit(S.fitted_f,S.fitted_alpha,s),[1,1]);
S.f_list = linspace(0.01,10,100);
S.yth = powerfun(f_list,l1(1,:)); % fitted exponential function

l1(2,:) = fminsearch(@(s)powerfit(S2.fitted_f,S2.fitted_alpha,s),[1,1]);

c = {[0.8500 0.3250 0.0980],[0.3010 0.7450 0.9330]};

attenuation_fig = figure;
loglog(S.fitted_f,S.fitted_alpha,'d','MarkerFaceColor',c{1},'MarkerEdgeColor','black');
hold on
plot(S.f_list,S.yth,'--','Color',c{1},'LineWidth',1.5);
hold on 
loglog(S2.fitted_f,S2.fitted_alpha,'o','MarkerFaceColor',c{2},'MarkerEdgeColor','black');
hold on 
plot(S2.f_list,S2.yth,'--','Color',c{2},'LineWidth',1.5);
xlabel('$f \: \rm (Hz)$','Interpreter','latex');
ylabel('$\alpha \: \rm (m^{-1})$','Interpreter','latex');
grid on 
axis([4e-2 4 1e-3 1])
ax = gca;
ax.FontSize = 13;

power_law_txt(1,:) = ['$\alpha(f) = ' sprintf('%0.2f',l1(1,2)) 'f^{' sprintf('%0.2f',l1(1,1)) '}$'];
power_law_txt(2,:) = ['$\alpha(f) = ' sprintf('%0.2f',l1(2,2)) 'f^{' sprintf('%0.2f',l1(2,1)) '}$'];
% legend('Data','Fitted Data',power_law_txt,'Interpreter','latex','Location','southeast','FontSize',13)
legend('',power_law_txt(1,:),'',power_law_txt(2,:),'Interpreter','latex','Location','northwest','FontSize',15)
set_Papermode(gcf);

figname = 'Comparison_power_laws_0226_waves005_03102023_DJI0402';
fig_folder = 'C:/Users/sebas/OneDrive/Bureau/These PMMH/Conferences/MePhy_2024/Figures/';
saveas(gcf,[fig_folder figname],'fig')
saveas(gcf,[fig_folder figname],'png')
saveas(gcf,[fig_folder figname],'pdf')
%% FUNCTION SECTION 

function yth = powerfun(x,l)
    yth = l(2)*x.^l(1);
end 

function d = powerfit(x,y,l)
    yth = powerfun(x,l);
    d = sum((yth-y).^2);
end 
