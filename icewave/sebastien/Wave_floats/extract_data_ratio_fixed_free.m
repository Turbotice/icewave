
%% Extract data from curves on the graph

ax = gca; 
h = findobj(gca,'Type','line');

curve_1 = struct();
curve_1.x = h(2).XData;
curve_1.y = h(2).YData;

curve_2 = struct();
curve_2.x = h(3).XData;
curve_2.y = h(3).YData;

fixed = curve_1;
free = curve_2;

data = struct();
data.fixed = fixed;
data.free = free;

base = 'C:/Users/sebas/OneDrive/Bureau/These PMMH/Waves_float/Data/';
filename = [base 'ratios_Ao_Aw_fixed_free'];

save(filename,'data','-v7.3')


%% Check which curve corresponds to what

figure,
plot(curve_1.x,curve_1.y)
hold on 
plot(curve_2.x,curve_2.y)

figure,
plot(curve_1.x,curve_1.y./curve_2.y)
