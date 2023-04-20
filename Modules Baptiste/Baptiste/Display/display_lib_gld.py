import matplotlib.pyplot as plt 
import matplotlib as mpl
import numpy as np

## TAILLE DES FIGURES

def set_size(width, fraction=1, subplots=(1, 1)):
    """Set figure dimensions to avoid scaling in LaTeX.

    Parameters
    ----------
    width: float or string
            Document width in points, or string of predined document type
    fraction: float, optional
            Fraction of the width which you wish the figure to occupy
    subplots: array-like, optional
            The number of rows and columns of subplots.
    Returns
    -------
    fig_dim: tuple
            Dimensions of figure in inches
    """
    if width == 'thesis':
        width_pt = 426.79135
    elif width == 'beamer':
        width_pt = 307.28987
    else:
        width_pt = width

    # Width of figure (in pts)
    fig_width_pt = width_pt * fraction
    # Convert from pt to inches
    inches_per_pt = 1 / 72.27

    # Golden ratio to set aesthetic figure height
    # https://disq.us/p/2940ij3
    #golden_ratio = (5**.5 - 1) / 2
    golden_ratio = .9
    # Figure width in inches
    fig_width_in = fig_width_pt * inches_per_pt
    # Figure height in inches
    fig_height_in = fig_width_in * golden_ratio * (subplots[0] / subplots[1])
    #fig_height_in = fig_width_in * .9 * (subplots[0] / subplots[1])

    return (fig_width_in, fig_height_in)


## FONTS

tex_fonts = {
    # Use LaTeX to write all text
    "text.usetex": True,
    "font.family": "serif",
    # Use 10pt font in plots, to match 10pt font in document
    "axes.labelsize": 14,
    "font.size": 14,
    # Make the legend/label fonts a little smaller
    "legend.fontsize":11,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12
}

plt.rcParams.update(tex_fonts)
mpl.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble'] = [
r'\usepackage{lmodern}', #lmodern: lateX font; tgheros: helvetica font; helvet pour helvetica
r'\usepackage{sansmath}', # math-font matching helvetica
r'\sansmath' # actually tell tex to use it!
r'\usepackage[scientific-notation=false]{siunitx}', # micro symbols
r'\sisetup{detect-all}', # force siunitx to use the fonts
]

## CHOIX DES COULEURS
n = 10
vcolors = plt.cm.viridis(np.linspace(0,1,n))

# def figurejolie(subplot = False, posplot = [1,1]):
#     if subplot == False :
#         plt.figure(figsize = set_size(width=350, fraction = 1, subplots = (1,1)))
        
#     else :
#         fig = plt.subplots(subplot[0],subplot[1], figsize = set_size(width = 700, subplots=subplot))
#         axes = []
                           
#         return fig, axes
    
def figurejolie(subplot = False, num_fig = None):
    if subplot == False :
        plt.figure(num = num_fig, figsize = set_size(width=500, fraction = 1, subplots = (1,1)))
        
    else :
        fig = plt.figure(num = num_fig, figsize = set_size(width=1000, fraction = 1, subplots = subplot))
        axes = []
                           
        return fig, axes


def joliplot(xlabel, ylabel, xdata, ydata, color = False, fig = False, axes = [], title = False, subplot = False, legend = False, log = False, exp = True, image = False, zeros = False):
    
    """Jolis graphs predefinis"""
    n = 16
    markers = ['0','x','o','v','p','X','d','s','s','h','.','.','o','o','o','v','v']
    markeredgewidth = [1.5,1.8,1.5, 1.5, 2.5, 1.3, 1.3, 1.6,1.6,2,1.6,1.6,2,2,2,2,2]
    ms = [7,6.5,7, 7, 9.2, 8, 8, 8, 7,7, 7, 7,9,7,7,7,7]
    mfc = ['None','#91A052','#990000',vcolors[5],'None','None','None','None','k','None','None','None','None','None','None','None','None']
    colors = ['g','#91A052','#990000', vcolors[5], '#008B8B', vcolors[2], '#FF8000', vcolors[6], 'k',vcolors[1],'#01FA22',vcolors[3], vcolors[1],'#990000', vcolors[2],'#990000', vcolors[2] ]
    
    """Pour un simple plot"""
    if subplot == False:
        
        #si c'est une image la mettre dans image
        if type(image) != bool :
            plt.imshow(image, cmap = plt.cm.gray)
            plt.xlabel(xlabel) 
            plt.ylabel(ylabel)
            if title != False:
                plt.title(title)
            plt.grid('off')
            plt.axis('off')
            plt.tight_layout()
        
        #pour un graph
        else :
            #si color = False fait une couleur au hasard
            if color == False:
                color = np.random.randint(1,n+1)
            
                
            if exp :
                marker = ' ' + markers[color]
            
            else :
                marker = '-'
    
            
            
            if title != False:
                plt.title(title)
               
            if legend != False :
                plt.plot(xdata, ydata, marker, color = colors[color], mfc = mfc[color], markeredgewidth = markeredgewidth[color], ms = ms[color], label = legend)
                plt.legend()
            else :
                plt.plot(xdata, ydata, marker, color = colors[color], mfc = mfc[color], markeredgewidth = markeredgewidth[color], ms = ms[color])
                
            plt.xlabel(xlabel) 
            plt.ylabel(ylabel)
            if log :
                plt.yscale('log')
                plt.xscale('log')
                plt.tight_layout()
                
            plt.grid()
            
            if zeros :
               plt.xlim(left=0 )
               plt.ylim(bottom=0)
                
    
    #pour un subplot
    else:
        if image != False :
            axes.append( fig.add_subplot(subplot[0], subplot[1], len(axes) + 1) )
            axes[-1].imshow(image, cmap = plt.cm.gray)
            axes[-1].set_xlabel(xlabel)
            axes[-1].set_ylabel(ylabel)
            if title != False:
                axes[-1].set_title(title)
            axes[-1].grid('off')
            plt.tight_layout()
                
        else :
            if color == False:
                color = np.random.randint(1,n+1)
                    
            if exp :
                marker = ' ' + markers[color]
            if exp == False :
                marker = '-'
    
                
            axes.append( fig.add_subplot(subplot[0], subplot[1], len(axes) + 1) )
            axes[-1].plot(xdata, ydata, marker, color = colors[color], mfc = mfc[color], markeredgewidth = markeredgewidth[color], ms = ms[color], label = legend)
            axes[-1].set_xlabel(xlabel)
            axes[-1].set_ylabel(ylabel) 
                
            if title != False :
                axes[-1].set_title(title)
            if legend != False :
                axes[-1].legend()
            if log == True :
                axes[-1].set_yscale('log')
                axes[-1].set_xscale('log')
                plt.tight_layout()
            
        return axes
            
                









