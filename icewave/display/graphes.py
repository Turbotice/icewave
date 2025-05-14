

import matplotlib.pyplot as plt
import numpy as np
import os
import pickle



def legende(x_legend,y_legend,title,ax=None,font=18,display=False,cplot=False,show=True,**kwargs):
    """
    Add a legend to the current figure
        Contains standard used font and sizes
        return a default name for the figure, based on x and y legends
    INPUT
    -----
    x_legend : str
        x label
    y_legend : str
        y label
    title : str
        title label
    colorplot : bool
        default False
        if True, use the title for the output legend
    OUTPUT
    -----
    fig : dict
        one element dictionnary with key the current figure number
        contain a default filename generated from the labels
    """
    if ax is None:
        ax = plt.gca()
    else:
        plt.sca(ax)
    #additionnal options ?
    plt.rc('font',family='Times New Roman')
    plt.rc('text', usetex=False)
    
    if r'\frac' in x_legend:
        plt.xlabel(x_legend,fontsize=font+4)
    else:
        plt.xlabel(x_legend,fontsize=font)
        
    if r'\frac' in y_legend:
        plt.ylabel(y_legend,fontsize=font+4,rotation = 0)
        ax.yaxis.set_label_coords(-0.13, 0.5) 
    else:
        plt.ylabel(y_legend,fontsize=font,*kwargs)#,rotation = 0)        
#    plt.ylabel(y_legend,fontsize=font)
    plt.title(title,fontsize=font)
    
    if show:
        refresh()

    #rec
    fig = figure_label(x_legend,y_legend,title,display=display,cplot=cplot) #fig is a dictionnary where the key correspond to the fig number and the element to the automatic title
    fig = get_data(fig)

    return fig


def make_title(keylist,d):
    s = r''
    for key in keylist:
        val = str(np.round(d[key]['value'],decimals=1))
        s = s+r'$'+d[key]['legend']+'$ = '+val+d[key]['unit']+', '
    s = s[:-2]
    return s

def get_data(fig,cplot=False):
    current = plt.gcf()
    lines = plt.gca().get_lines()
    
    Dict = {}
    for i,line in enumerate(lines):
        xd = line.get_xdata()
        yd = line.get_ydata()
        Dict['xdata_'+str(i)]=xd
        Dict['ydata_'+str(i)]=yd
        
        if cplot:
            zd = line.get_zdata()
            Dict['zdata'+str(i)]=zd
                
    fig[current.number]['data']=Dict  
    return fig
    
def figure_label(x_legend,y_legend,title,display=True,cplot=False,include_title=False):    
    #generate a standard name based on x and y legend, to be used by default as a file name output
    x_legend = remove_special_chars(x_legend)
    y_legend = remove_special_chars(y_legend)
    title = remove_special_chars(title)
    
    fig={}
    current = plt.gcf()
    fig[current.number]={}
    if cplot:
        fig[current.number]['fignum']=title#+'_'+x_legend+'_'+y_legend #start from the plotted variable (y axis)
    else:
        fig[current.number]['fignum']=y_legend+'_vs_'+x_legend #start from the plotted variable (y axis)
        
    if include_title:
        fig[current.number]['fignum']=y_legend+'_vs_'+x_legend+'_'+title
        
    if display:
        print(current.number,fig[current.number])
    return fig


def plot(fun,x,y,fignum=1,label='-',subplot=None,**kwargs):
    """
    plot a graph using the function fun
    fignum can be specified
    any kwargs from plot can be passed
    Use the homemade function refresh() to draw and plot the figure, no matter the way python is called (terminal, script, notebook)
    """
    #set_fig(fignum,subplot=subplot)
    fun(x,y,label,**kwargs)
    refresh()   

def graph(x,y,fignum=1,label='-',subplot=None,**kwargs):
    """
    plot a graph using matplotlib.pyplot.plot function
    fignum can be specified
    cut x data if longer than y data
    any kwargs from plot can be passed
    Use the homemade function refresh() to draw and plot the figure, no matter the way python is called (terminal, script, notebook)
    """    
    xp = np.asarray(x)
    yp = np.asarray(y)
    if len(xp)>len(yp):
        print("Warning : x and y data do not have the same length")
        xp=xp[:len(yp)]    
        
    plot(plt.plot,xp,yp,fignum=fignum,label=label,subplot=subplot,**kwargs)
    
def errorbar(x,y,xerr,yerr,fignum=1,label='k^',subplot=None,**kwargs):
    """
    plot a graph using matplotlib.pyplot.errorbar function
    fignum can be specified
    cut x data if longer than y data
    any kwargs from plot can be passed
    """
    set_fig(fignum,subplot=subplot)        
    plt.errorbar(x,y,yerr,xerr,label,**kwargs)
    refresh()

def set_fig(fignum,subplot=None):
    if fignum==-2:
        fig=None
        pass
    if fignum==-1:
        fig=plt.figure()
    if fignum==0:
        fig=plt.cla()
    if fignum>0:
        fig=plt.figure(fignum)
        
    if subplot is not None:
        # a triplet is expected !
        ax = fig.add_subplot(subplot)
        return fig,ax
    return fig
    
def cla(fignum):
    set_fig(fignum)
    plt.cla()
    
def set_axis(xmin,xmax,ymin,ymax):
    #set axis
    plt.xlim(xmin=xmin,xmax=xmax)
    plt.ylim(ymin=ymin,ymax=ymax)        
    
def time_label(M,frame):
    Dt = M.t[frame+1]-M.t[frame]
    title = 't = '+str(floor(M.t[frame]*1000)/1000.)+' s, Dt = '+str(floor(Dt*10**4)/10.)+' ms'
    return title
    
def pause(time=3):
    plt.pause(time)
    
def refresh(hold=True,block=False,ipython=True):
    """
    Refresh the display of the current figure. 

    INPUT 
    hold (opt): False if the display has overwritten.

    OUTPUT 
    None
    """    
    if not ipython:           
#        plt.pause(0.1)
#        plt.draw()   
        plt.hold(hold)
        plt.show(block)
 
def remove_special_chars(string,chars_rm=['$','\ ','[',']','^','/',') ','} ',' '],chars_rp=['{','(',',','=','.']):
    """
    Remove characters from a typical latex format to match non special character standards
    INPUT
    -----
    string : str
        input string
    chars_rm : str list. Default value : ['$','\ ',') ']
        char list to be simply removed
    chars_rp : str list. Default value : ['( ',',']
        char list to be replaced by a '_'
    OUTPUT
    -----
    string : str
        modified string
    """
    for char in chars_rm:
        string = string.replace(char[0],'')
    for char in chars_rp:
        string = string.replace(char[0],'_')
    return string
    
def colorbar(fignum=-2,label=''):
    set_fig(fignum)
    c=plt.colorbar()        
    c.set_label(label,fontsize=18)
    return c
    
def set_title(M,opt=''):
    #if Zplane attribute exist !
    title='Z= '+str(int(M.param.Zplane))+', '+M.param.typeview+', mm, '+M.Id.get_id()+', '+opt 
    plt.title(title,fontsize=18)
    return title

def clegende(c,c_legend):
    c.set_label(c_legend)    

def save_graphes(M,figs,prefix='',suffix=''):
    save_figs(figs,savedir='./Results/'+os.path.basename(M.dataDir)+'/',prefix=prefix,suffix=suffix)

def save_figs(figs,savedir='./',suffix='',prefix='',frmt='pdf',dpi=300,display=False,overwrite=True):
    """
    save a dictionnary of labeled figures using dictionnary elements 
        dict can be autogenerated from the output of the graphes.legende() function
        by default the figures are saved whithin the same folder from where the python code has been called
    INPUT
    -----
    figs : dict of shape {int:str,...}
        the keys correspond to the figure numbers, the associated field
    savedir : str. default : './'
    frmt : str
        file Format
    dpi : int
        division per inch. set the quality of the output
    OUTPUT
    -----
    None
    """
    c=0
    filename=''
    for key in figs.keys():
        fig = figs[key]        
        #save the figure
        filename = savedir+prefix+fig['fignum']+suffix
        save_fig(key,filename,frmt=frmt,dpi=dpi,overwrite=overwrite)
        c+=1

        #save the data
        if not fig['data']=={}:
            with open(filename+'.pkl', 'wb') as handle:
                pickle.dump(fig['data'], handle, protocol=pickle.HIGHEST_PROTOCOL)
            #print(fig['data'])
            #try:
            #    h5py_s.save(filename,fig['data'])
            #except:
            #    print('Data not saved')
    if display:
        print('Number of auto-saved graphs : '+str(c))
        print(filename)
    
def save_fig(fignum,filename,frmt='pdf',dpi=300,overwrite=False):
    """
    Save the figure fignumber in the given filename
    INPUT
    -----
    fignumber : number of the fig to save
    filename : name of the saving file
    fileFormat (opt)
    dpi (opt) : number of dpi for image based format. Default is 300
    OUTPUT 
    -----
    None
    """
    if os.path.dirname(filename)[0]=='.':
        Dir = '/Users/stephane/Documents/git/stephane/stephane/iPython_notebooks'+os.path.dirname(filename)[1:]
        filename = Dir +'/'+ os.path.basename(filename)
    else:
        Dir = os.path.dirname(filename)
    if not os.path.isdir(Dir):
        os.makedirs(Dir)

    filename=filename+'.'+frmt
    if fignum!=0:
        plt.figure(fignum)    
    
    if not os.path.isfile(filename):
        plt.savefig(filename, dpi=dpi)
    elif overwrite:
        plt.savefig(filename, dpi=dpi)
    else:
        print("figure already exists")
        
def plot_axes(fig,num):
    ax = fig.add_subplot(num)
    ax.set_aspect('equal', adjustable='box')
    #graphes.legende('','','Front view')
    #draw_fieldofview(M.Sdata,ax3,view='front')
    return ax
    
def color_plot(x,y,Z,fignum=1,vmin=0,vmax=0,log=False,show=False,cbar=False):
    """
    Color coded plot
    INPUT
    -----	
    x : 2d numpy array
    y : 2d numpy array
    Z : 2d numpy array 
    OUTPUT
    ------
    None
    	"""
    fig,ax = set_fig(fignum,subplot=111)
    #ax.axis('off')
        
    if log:
        Z = np.log10(Z)
         
    if vmin==vmax==0:
        c=plt.pcolormesh(x,y,Z)#,shading='gouraud')
    else:
        c=plt.pcolormesh(x,y,Z,vmin=vmin,vmax=vmax)

    if cbar:
        colorbar()
    if show:
        refresh()
  #  return fig,ax,c
    
def get_axis_coord(M,direction='v'):
    X = M.x
    Y = M.y
    
    if hasattr(M,'param'):
        if M.param.angle==90:
            Xb = X
            Yb = Y
            X = Yb
            Y = Xb
        
    return X,Y    
#    return rotate(X,Y,M.param.angle)
    
def rotate(X,Y,angle):
    angle = angle/180*np.pi
    return X*np.cos(angle)-Y*np.sin(angle),Y*np.cos(angle)+X*np.sin(angle)


######################################################################
#################### Histograms and pdfs #############################
######################################################################

def hist(Y,Nvec=1,fignum=1,num=100,step=None,label='o-',log=False,normalize=True,xfactor=1,**kwargs):
    """
    Plot histogramm of Y values
    """
    set_fig(fignum)
   # print('Number of elements :'+str(len(Y)))
    if step is None:
        n,bins=np.histogram(np.asarray(Y),bins=num,**kwargs)
      #  print(bins)
    else:
        d = len(np.shape(Y))
#        print('Dimension : '+str(d))
        N=np.prod(np.shape(Y))
        if N<step:
            step=N
        n,bins=np.histogram(np.asarray(Y),bins=int(N/step))
        
    if normalize:
        dx = np.mean(np.diff(bins))
        n=n/(np.sum(n)*dx)
        
    xbin=(bins[0:-1]+bins[1:])/2/xfactor
    n = n*xfactor
       
    if log:
    # Plot in semilogy plot
        semilogy(xbin/Nvec,n,fignum,label)
    else:
        plt.plot(xbin,n,label)
        plt.axis([np.min(xbin),np.max(xbin),0,np.max(n)*1.1])
        
    refresh()
    return xbin,n    
    
def pdf(M,field,frame,Dt=10,Dx=1024,label='ko-',fignum=1,a=15.,norm=True,sign=1):
    import stephane.manager.access as access
    Up = access.get(M,field,frame,Dt=Dt)
    
    limits = [(0,Dx),(0,Dx)]
    Up = sign*access.get_cut(M,field,limits,frame,Dt=Dt)
    
    figs = distribution(Up,normfactor=1,a=a,label=label,fignum=fignum,norm=norm)
    
    return figs

def pdf_ensemble(Mlist,field,frame,Dt=10,Dx=1024,label='r-',fignum=1,a=10.,norm=True,model=False):

    import stephane.manager.access as access

    U_tot = []
    
    for M in Mlist:
        pdf(M,field,frame,Dt=Dt,Dx=Dx,label='k',fignum=fignum,a=a,norm=False)
        
        Up = access.get(M,field,frame,Dt=Dt)
       # limits = [(0,Dx),(0,Dx)]
    #    Up = access.get_cut(M,field,limits,frame,Dt=Dt) 
        # if Dx is larger than the box size, just keep all the data
        U_tot = U_tot + np.ndarray.tolist(Up)
        
    N = len(Mlist)
    U_tot = np.asarray(U_tot)
    
    x,y,figs = distribution(U_tot,normfactor=N,a=a,label=label,fignum=fignum,norm=norm)
    
    if model:
        n = len(y)
        b = y[n//2]    
        Dy = np.log((y[n//2+n//8] + y[n//2-n//8])/2./b)
    
        a = - Dy/x[n//2+n//8]**2
    
        P = b*np.exp(-a*x**2)
        semilogy(x,P,label='b.-',fignum=fignum)
    
    set_axis(min(x),max(x),1,max(y)*2)
    if field=='omega' or field=='strain':
        unit = ' (s^-1)'
    elif field=='E':
        unit = 'mm^2/s^2'
    else:
        unit = ' (mm/s)'
    figs = {}
    figs.update(legende(field+unit,field+' PDF',time_label(M,frame)))
    return figs
    
def avg_from_dict(dd,keyx,keyy,times,fignum=1,display=True,label='b-'):
    """
    Compute the average function from a dictionnary with keys (time,keyx) (time,keyy) for time in times
    
    """    
    avg = {}
    avg[keyx]= np.mean([dd[(time,keyx)] for time in times],axis=0)
    avg[keyy]= np.mean([dd[(time,keyy)] for time in times],axis=0)
    
    std = {}
    std[keyx]= np.std([dd[(time,keyx)] for time in times],axis=0)
    std[keyy]= np.std([dd[(time,keyy)] for time in times],axis=0)
    
    if display:
        for time in times:
            graph(dd[(time,keyx)],dd[(time,keyy)],label='k-',fignum=fignum,color='0.7')
                        
        errorbar(avg[keyx],avg[keyy],std[keyx],std[keyy],fignum=fignum,label=label)
    
    return avg,std

def distribution(Y,normfactor=1,a=10.,label='k',fignum=1,norm=True):

    Y = np.asarray(Y)
    Median = np.sqrt(np.nanmedian(Y**2))
    
    "test if the field is positive definite"
    t = Y>=0
    
    if norm:
        bound = a
    else:
        bound = a*Median
    step = bound/10**2.5

    if t.all():
        x = np.arange(0,bound,step)
    else:
        x = np.arange(-bound,bound,step)
    
    if norm:
        Y = Y/Median
        
    n,bins = np.histogram(Y,bins=x)
    xbin = (bins[:-1]+bins[1:])/2
    n = n/normfactor # in case of several series (ensemble average)

    semilogy(xbin,n,label=label,fignum=fignum)
    set_axis(min(xbin),max(xbin),min(n)/2,max(n)*2)
    figs = {}
    figs.update(legende('','PDF',''))
    
    val = 0.5
    x_center = xbin[np.abs(xbin)<val]
    n_center = n[np.abs(xbin)<val]
    
    moy = np.sum(n*xbin)/np.sum(n)
    std = np.sum(n*(xbin-moy)**2)/np.sum(n)
    
    print("Distribution : "+str(moy)+' +/- '+str(std))
    #    a = fitting.fit(fitting.parabola,x_center,n_center)
    #    n_th = fitting.parabola(xbin,a)
    #    graph(xbin,n_th,label='r-',fignum=fignum)
    return xbin,n,figs
