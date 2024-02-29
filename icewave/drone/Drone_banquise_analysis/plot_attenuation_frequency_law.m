
close all;
clear all;

%% Load data
base = '//192.168.1.70/Share/Data/0215/Drones/mesange/matData/waves_003/';
filename = 'attenuation_coef_data_fmin0p07_fmax1p1_f_thresh0p8.mat';

fullname = [base filename];

fig_folder = [base 'Plots/'];

%%
disp('Loading data...')
load(fullname)
disp('Data loaded')

%%

mask = dist_fit < 0.2 ;
f = freq';
fitted_f = f(mask);
fitted_alpha = lambda(mask);

save_boolean = 0;

l1 = fminsearch(@(s)powerfit(fitted_f,fitted_alpha,s),[1,1]);
% l1 = fminsearch(@(s)powerfit(fitted_k,fitted_omega,s),[1,1]);
f_list = linspace(0.01,2,100);
yth = powerfun(f_list,l1); % fitted exponential function

% Take more datas 

% theory = powerfun(freq_secondary,l1);
% Need to take the orthogonal distance in the logarithm scale....
% dist_to_fit = sum((log((theory)) - log((lambda_secondary))).^2)/sum(log(abs(lambda_secondary)).^2);
% 
% mask = dist_to_fit < 0.5;
% freq_plot = freq_secondary(mask);
% lambda_plot = lambda_secondary(mask);

attenuation_fig = figure;
% loglog(freq_secondary,lambda_secondary,'o');
% hold on 
loglog(fitted_f,fitted_alpha,'o','MarkerFaceColor',[0.3010 0.7450 0.9330],'MarkerEdgeColor','black');
hold on
plot(f_list,yth,'r--','LineWidth',1.5);
% plot(freq_secondary,theory,'k-')
xlabel('$f \: \rm (Hz)$','Interpreter','latex');
ylabel('$\alpha \: \rm (m^{-1})$','Interpreter','latex');
grid on 
axis([0.05 2 0.0001 1])
ax = gca;
ax.FontSize = 13;

power_law_txt = ['$\alpha(f) = ' sprintf('%0.2f',l1(2)) 'f^{' sprintf('%0.2f',l1(1)) '}$'];
% legend('Data','Fitted Data',power_law_txt,'Interpreter','latex','Location','southeast','FontSize',13)
legend('Data',power_law_txt,'Interpreter','latex','Location','southeast','FontSize',13)

set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'filename','-dpdf','-r0')

if save_boolean
    attenuation_filename = [fig_folder 'attenuation_law'];
    saveas(attenuation_fig,attenuation_filename,'fig');
    saveas(attenuation_fig,attenuation_filename,'pdf');
end 

%%

function yth = powerfun(x,l)
    yth = l(2)*x.^l(1);
end 

function d = powerfit(x,y,l)
    yth = powerfun(x,l);
    d = sum((yth-y).^2);
end 