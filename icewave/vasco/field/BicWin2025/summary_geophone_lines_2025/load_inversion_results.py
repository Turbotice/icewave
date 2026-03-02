#%%
import pickle
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import csv

#%%

def open_inverted_data(pkl_file_path):

    with open(pkl_file_path,'rb') as file:
        data = pickle.load(file)


    T_critic = data[0]
    X = data[1]
    misfit_accepted = data[2]

    nonzeroidx = np.nonzero(X[0, :])[0][-1]

    data_to_plot = [X[i, :nonzeroidx] for i in range(4)]

    # Step 3: Create a 2x2 subplot and plot histograms with KDE
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    titles = ['Thickness (m)', 'E (Pa)',
            'nu ', 'rho (kg/m3)']

    dd = {}

    for i, ax in enumerate(axes.flat):
        data = data_to_plot[i]
        
        dd['data_to_plot_'+str(i)] = {'values':data}

        # Plot the histogram
        sns.histplot(data, bins=30, kde=True, ax=ax, edgecolor='black', color='blue', alpha=0.6)

        # Set titles and labels
        ax.set_title(titles[i])
        ax.set_xlabel('Value')
        ax.set_ylabel('PDF')

    plt.tight_layout()
    plt.show()

    return dd

#%% Try to find maximum of kde
from scipy import stats

def find_max_kde(data_to_plot_i):

    kde = stats.gaussian_kde(data_to_plot_i)
    Y = kde.evaluate(data_to_plot_i)

    # compute argmax
    return np.argmax(Y) # this will index the data_to_plot_i array (i.e. same indices than the original data)



# %% postprocess all
summary_dir = "B:/General/Summary_geophone_lines/"
txtfile_inverted_acq = summary_dir + "inversions/acquisitions_nb_inverted.txt"

dict_results = {}

dict_results['inverted_data'] = np.loadtxt(txtfile_inverted_acq,dtype=str)

# %%
inverted_acqs = []
inverted_dates = []
for i in range(len(dict_results['inverted_data'])):
    inverted_acqs.append(int(dict_results['inverted_data'][i][-1]))
    inverted_dates.append(dict_results['inverted_data'][i][:4])

print(inverted_acqs)
print(inverted_dates)
# %%
disk = 'B:'
year = 2025

i=0
idx_datatoplot = 1 # thickness : 0 ; Young: 1 ; Poisson : 2 ; rho : 3
dateacq = dict_results['inverted_data'][i]
date = inverted_dates[i]
acqstr = str(inverted_acqs[i]).zfill(4)

pkl_file_path = f"{disk}/Data/{date}/Geophones/{year}_{date}_acq{acqstr}_bidir_inversion.pkl"

# test :

dict_results[dateacq] = open_inverted_data(pkl_file_path)

argmaxkde = find_max_kde(dict_results[dateacq]['data_to_plot_'+str(idx_datatoplot)]['values'])
maxvalkde = dict_results[dateacq]['data_to_plot_'+str(idx_datatoplot)]['values'][argmaxkde]



# %% range tous les histogrammes dans un dico

for i in range(len(dict_results['inverted_data'])):
    dateacq = dict_results['inverted_data'][i]
    date = inverted_dates[i]
    acqstr = str(inverted_acqs[i]).zfill(4)

    pkl_file_path = f"{disk}/Data/{date}/Geophones/{year}_{date}_acq{acqstr}_bidir_inversion.pkl"

    dict_results[dateacq] = open_inverted_data(pkl_file_path)

#%% ranger les max des kde por chaque quantité (epaisseur, young...)

list_argmaxkde = []
list_maxvalkde =[]

argmaxkde = find_max_kde(dict_results[dateacq]['data_to_plot_'+str(idx_datatoplot)]['values'])
maxvalkde = dict_results[dateacq]['data_to_plot_'+str(idx_datatoplot)]['values'][argmaxkde]

# trouver le maxkde (meilleure estimation) pour chaque param de glace

dict_results[dateacq]

for i in range(len(dict_results['inverted_data'])):
    dateacq = dict_results['inverted_data'][i]
    date = inverted_dates[i]
    for idx_datatoplot in [0,1,2,3]:
        argmaxkde = find_max_kde(dict_results[dateacq]['data_to_plot_'+str(idx_datatoplot)]['values'])
        maxvalkde = dict_results[dateacq]['data_to_plot_'+str(idx_datatoplot)]['values'][argmaxkde]
        dict_results[dateacq]['data_to_plot_'+str(idx_datatoplot)]['argmaxkde'] = argmaxkde
        dict_results[dateacq]['data_to_plot_'+str(idx_datatoplot)]['maxvalkde'] = maxvalkde

# %% create csv file with all information (including coordinates and water depths)
coord_water_level_filepath = f"{disk}/General/Summary_geophone_lines/water_level_acquisitions/coord_waterlevel_allacq.csv"

with open(coord_water_level_filepath, mode ='r')as file:
    csvFile = csv.reader(file, delimiter=',')
    count=0
    data = []
    for lines in csvFile:
        if count==0:
            header = lines
        else:
            data.append(lines)
        count+=1
print(header)

data = np.array(data)


data_new = np.empty((data.shape[0], data.shape[1]+4),dtype=object)
data_new[:data.shape[0],:data.shape[1]] = data
header_new = header + ['thickness (m)', 'Young modulus (Pa)', 'Poisson ratio', 'density (kg/m^3)']

for i in range(len(dict_results['inverted_data'])):
    dateacq = dict_results['inverted_data'][i]
    date = inverted_dates[i]
    completedatestr = f'{year}{date}'
    acqnum_str = str(inverted_acqs[i])
    idx_row = np.where((data[:,0]==completedatestr)&(data[:,1]==acqnum_str))[0][0]
    print(idx_row)
    for idx_datatoplot in [0,1,2,3]:
        data_new[idx_row,data.shape[1]+idx_datatoplot] = dict_results[dateacq]['data_to_plot_'+str(idx_datatoplot)]['maxvalkde']

csvpath2save = f"{disk}/General/Summary_geophone_lines/results_acquisitions.csv"

with open(csvpath2save, 'w', encoding='UTF8', newline='') as f:
    writer = csv.writer(f,delimiter=',')

    # write the header
    writer.writerow(header_new)

    # write multiple rows
    writer.writerows(data_new)
# %%
