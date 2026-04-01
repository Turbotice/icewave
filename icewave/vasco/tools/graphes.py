#%%
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import os
import re
import pickle
import matplotlib.lines as mlines
#%%
def defaults_params_errorbar_colorvar(dir2save='C:/Users/Vasco Zanchi/Desktop/AutoOutputFigures'):
    dict_params = {}
    dict_params['xerr'] = None
    dict_params['yerr'] = None
    dict_params['xlabel'] = None
    dict_params['ylabel'] = None
    dict_params['z'] = None
    dict_params['zlabel'] = None
    dict_params['vmin'] = None
    dict_params['vmax'] = None
    dict_params['symbols'] = 'o'
    dict_params['title'] = None
    dict_params['savefig'] = False
    dict_params['savedata'] = False
    dict_params['dir2save'] = dir2save
    dict_params['capsize'] = 2
    dict_params['ecolor'] = 'gray'
    dict_params['xscale'] = 'linear'
    dict_params['yscale'] = 'linear'
    dict_params['xlim'] = None
    dict_params['ylim'] = None
    dict_params['fontsize'] = 15
    dict_params['cmap'] = 'viridis'
    dict_params['alpha'] = 1.0
    dict_params['figsize'] = (6, 4)
    
    return dict_params



def errorbar_propre_colorvar(dict_params=dict):
    """
    Ici on trace un errorbar, avec comme couleurs de 
    points une 3eme variable z -> variation continue 
    des couleurs sur une échelle entre vmin et vmax
    """
    x = dict_params['x']
    y = dict_params['y']
    xerr = dict_params['xerr']
    yerr = dict_params['yerr']
    xlabel = dict_params['xlabel']
    ylabel = dict_params['ylabel']
    z = dict_params['z'] # 3eme variable, représentée en niveaux de couleurs
    zlabel = dict_params['zlabel']
    vmin = dict_params['vmin']
    vmax = dict_params['vmax']
    symbols = dict_params['symbols']
    title = dict_params['title']
    savefig = dict_params['savefig']
    savedata = dict_params['savedata']
    dir2save = dict_params['dir2save']
    capsize = dict_params['capsize']
    ecolor = dict_params['ecolor']
    xscale = dict_params['xscale']
    yscale = dict_params['yscale']
    xlim = dict_params['xlim']
    ylim = dict_params['ylim']
    fontsize = dict_params['fontsize']
    alpha = dict_params['alpha']
    cmap = dict_params['cmap']
    figsize =dict_params['figsize']

    today = datetime.now().strftime('%Y%m%d')
    timenow = datetime.now().strftime('%H%M%S')
    if (savedata|savefig)&(os.path.exists(f'{dir2save}/{today}')==False):
        os.mkdir(f'{dir2save}/{today}')
        print('created directory to save figure/data')
    figname = f"{timenow}_{ylabel or 'y'}_vs_{xlabel or 'x'}"
    figname = re.sub(r'[\\/*?:"<>|$\[\]\s]', '_', figname) # on enleve les caractere interdit pour la sauvegarde sur windows
    
    plt.figure(figsize=figsize)
    plt.errorbar(x, y,yerr=yerr,xerr=xerr,linestyle='',marker='',ecolor=ecolor, capsize=capsize)
    # Scatter
    sc = plt.scatter(
        x, y,
          c=z,
           cmap=cmap,
             marker=symbols,
               vmin=vmin,
                 vmax=vmax,
                   zorder=2
                   )

    # Récupérer le cmap et normalisation
    cmap = sc.cmap
    norm = sc.norm

    clb = plt.colorbar(sc)
    clb.ax.set_title(zlabel,fontsize=8)
    plt.xlabel(xlabel,fontsize=fontsize)
    plt.ylabel(ylabel,fontsize=fontsize)
    if xlim is not None:
        plt.xlim(xlim[0], xlim[1])
    if ylim is not None:
        plt.ylim(ylim[0], ylim[1])
    plt.yscale(yscale)
    plt.xscale(xscale)
    #plt.legend()
    if savefig:
        plt.savefig(f'{dir2save}/{today}/{figname}.pdf', dpi=300)
    if savedata:
        with open(f"{dir2save}/{today}/dict_params_{timenow}.pkl", 'wb') as f:
            pickle.dump(dict_params, f)

    plt.show()


# %% now functions to plot in the case where there are finite categories

def defaults_params_errorbar_categories(dir2save='C:/Users/Vasco Zanchi/Desktop/AutoOutputFigures'):
    dict_params = {}
    dict_params['xerr']      = None
    dict_params['yerr']      = None
    dict_params['xlabel']    = None
    dict_params['ylabel']    = None
    dict_params['z']         = None   # 1D array of category labels (e.g. ['A','A','B','C',...])
    dict_params['zlabel']    = None   # legend title
    dict_params['colors']    = None   # dict mapping category -> color, e.g. {'A':'red','B':'blue'}
    dict_params['symbols']   = None   # dict mapping category -> marker, e.g. {'A':'s','B':'o'}
    dict_params['title']     = None
    dict_params['savefig']   = False
    dict_params['savedata']  = False
    dict_params['dir2save']  = dir2save
    dict_params['capsize']   = 2
    dict_params['ecolor']    = 'gray'   # if None, each category uses its own color for errorbars
    dict_params['xscale']    = 'linear'
    dict_params['yscale']    = 'linear'
    dict_params['xlim']      = None
    dict_params['ylim']      = None
    dict_params['fontsize']  = 15
    dict_params['alpha']     = 1.0
    dict_params['figsize']   = (6, 4)
    dict_params['markersize']= 6
    return dict_params


def errorbar_categories(dict_params=dict):
    """
    Errorbar plot where points are colored and shaped by discrete category labels (z).
    Each unique category gets a unique color + symbol, shown in the legend.
    
    - z        : array-like of category labels, same length as x and y
    - colors   : dict {category: color}. If None, auto-assigned from tab10.
    - symbols  : dict {category: marker}. If None, cycles through a default list.
    - ecolor   : if None, errorbars take the color of their category.
    """


    x          = np.asarray(dict_params['x'])
    y          = np.asarray(dict_params['y'])
    xerr       = dict_params['xerr']
    yerr       = dict_params['yerr']
    xlabel     = dict_params['xlabel']
    ylabel     = dict_params['ylabel']
    z          = np.asarray(dict_params['z'])        # category labels
    zlabel     = dict_params['zlabel']               # legend title
    colors     = dict_params['colors']               # dict or None, example : {'a':'blue', 'b':'red'}
    symbols    = dict_params['symbols']              # dict or None
    title      = dict_params['title']
    savefig    = dict_params['savefig']
    savedata   = dict_params['savedata']
    dir2save   = dict_params['dir2save']
    capsize    = dict_params['capsize']
    ecolor     = dict_params['ecolor']
    xscale     = dict_params['xscale']
    yscale     = dict_params['yscale']
    xlim       = dict_params['xlim']
    ylim       = dict_params['ylim']
    fontsize   = dict_params['fontsize']
    alpha      = dict_params['alpha']
    figsize    = dict_params['figsize']
    markersize = dict_params['markersize']

    # --- Auto-assign colors and symbols if not provided ---
    categories = list(dict.fromkeys(z))  # unique values, order preserved
    default_markers = ['o', 's', '^', 'D', 'v', 'P', '*', 'X', 'h', '+']
    cmap_tab10 = plt.get_cmap('tab10')

    if colors is None:
        colors = {cat: cmap_tab10(i % 10) for i, cat in enumerate(categories)}
    if symbols is None:
        symbols = {cat: default_markers[i % len(default_markers)] for i, cat in enumerate(categories)}

    # --- Save setup ---
    today   = datetime.now().strftime('%Y%m%d')
    timenow = datetime.now().strftime('%H%M%S')
    if (savedata | savefig) and not os.path.exists(f'{dir2save}/{today}'):
        os.mkdir(f'{dir2save}/{today}')
        print('Created directory to save figure/data')
    figname = f"{timenow}_{ylabel or 'y'}_vs_{xlabel or 'x'}"
    figname = re.sub(r'[\\/*?:"<>|$\[\]\s]', '_', figname) # on enleve les caractere interdit pour la sauvegarde sur windows

    # --- Plot ---
    fig, ax = plt.subplots(figsize=figsize)

    for cat in categories:
        mask = (z == cat)
        xc   = x[mask]
        yc   = y[mask]
        xerrc = np.asarray(xerr)[mask] if xerr is not None else None
        yerrc = np.asarray(yerr)[mask] if yerr is not None else None
        color  = colors[cat]
        marker = symbols[cat]
        ec     = color if ecolor is None else ecolor  # errorbar color

        ax.errorbar(
            xc, yc,
            yerr=yerrc, xerr=xerrc,
            linestyle='', marker=marker,
            color=color, ecolor=ec,
            capsize=capsize, alpha=alpha,
            markersize=markersize,
            label=str(cat)
        )

    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.set_ylabel(ylabel, fontsize=fontsize)
    if title:
        ax.set_title(title, fontsize=fontsize)
    if xlim is not None:
        ax.set_xlim(xlim[0], xlim[1])
    if ylim is not None:
        ax.set_ylim(ylim[0], ylim[1])
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)

    legend = ax.legend(title=zlabel, fontsize=fontsize - 3, title_fontsize=fontsize - 3)

    plt.tight_layout()

    if savefig:
        plt.savefig(f'{dir2save}/{today}/{figname}.pdf', dpi=300)
    if savedata:
        with open(f"{dir2save}/{today}/dict_params_{timenow}.pkl", 'wb') as f:
            pickle.dump(dict_params, f)

    plt.show()
# %%
# exemple d'utilisation
"""
# ici pas de barres d'erreurs mais fonctionne aussi

dict_fig1 = graphes.defaults_params_errorbar_categories()

dict_fig1['x'] = arr_T_core_celsius
dict_fig1['y'] = arr_sigma_c
dict_fig1['z'] = np.empty(len(arr_T_core_celsius),dtype=object)
dict_fig1['z'][mask_HaHa] = 'HaHa'
dict_fig1['z'][mask_Hatee] = 'Hatee'
dict_fig1['z'][mask_Capelans_beach] = 'Capelans beach'
dict_fig1['z'][mask_Capelans_floating_ice] = 'Capelans floating'
dict_fig1['ylabel'] = '$\sigma_c$ [MPa]'
dict_fig1['xlabel'] = 'T [°C]'
dict_fig1['savefig'] = False
dict_fig1['title'] = 'sigma_c vs T for cores in different locations'
#dict_fig1['symbols'] = {'HaHa':'o','Hatee':'o','Capelans beach':'o','Capelans floating':'o'}


graphes.errorbar_categories(dict_params=dict_fig1)

"""