# ce code contient les fonctions utilisées dans la routine vitesse_phase
#  ainsi que la fonction routine vitesse phase elle-meme

####################################""""
# import modules
import numpy as np
import matplotlib.pyplot as plt
import h5py
import matplotlib
#import fitutils as fu
from scipy.signal import find_peaks
import os
import pickle
import csv
import re
from scipy.signal import savgol_filter
from scipy.optimize import curve_fit

# def fonctions
def savedict(filename,dico):
    with open(filename+'.pickle', 'wb') as handle:
        pickle.dump(dico, handle, protocol=pickle.HIGHEST_PROTOCOL)
def loaddict(filename):
    with open(filename+'.pickle', 'rb') as handle:
        dico = pickle.load(handle)
    return dico
def extract_data_mat_file(path_mat_file):
    with h5py.File(path_mat_file, 'r') as file:
        data = file['data'][:]  # Accessing 'data' variable
    return data

def load_complex_field(path_mat_file):
    matrix = extract_data_mat_file(path_mat_file)
    # Initialize an empty matrix to store complex numbers
    complex_matrix = np.empty(matrix.shape, dtype=complex)
    
    # Iterate through each element in the matrix
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            # Get the tuple from the original matrix
            a, b = matrix[i, j]
            
            # Create a complex number from the tuple and assign it to the corresponding position
            complex_matrix[i, j] = complex(a, b)
    return complex_matrix    

def load_csv_input_RDD(file_path,skipfirst2rows = True):
    # Dictionaries to store categories and lists
    list_categories = []
    dict_lists = {}
    list_types = []

    # Open and read the CSV file
    with open(file_path, 'r') as file:
        reader = csv.reader(file)
        count = 0

        for row in reader:
            if count == 0 and skipfirst2rows:  # Handle the header row
                list_categories = [row[i] for i in range(len(row))]  # Map index to column name
                for i in range(len(row)):
                    dict_lists[row[i]] = []  # Initialize empty lists for each column
                count += 1
            elif count == 1 and skipfirst2rows:
                list_types = [row[i] for i in range(len(row))]  # Map index to column name
                count+=1
            else:
                for i in range(len(row)):
                    if list_types[i]=='int':
                        dict_lists[list_categories[i]].append(int(row[i]))  # Append each cell to its corresponding list
                    elif list_types[i]=='float':
                        dict_lists[list_categories[i]].append(float(row[i]))
                    else: # par defaut ca sera str
                        dict_lists[list_categories[i]].append(row[i])
                count += 1
    
    return list_categories,dict_lists,list_types

##############################################

#  Routine permettant de calculer la vitesse de phase pour chaque freq

def routine_vitesse_phase(f_exc=float, freq_acq=float, general_folder=str, dosavgol=False,W=int,Dt=int,plot=True,direction=int,index_profile_line=int, xlim_fit = tuple, dcm = float, dpx = float):
    ## load data
    #general_folder = 'F:/Banquise/Vasco/Frigo_pmmh/' +date+ '/'
    #general_folder = f'F:/manip_grenoble2024/manips_relation_dispersion/{date}/Acquisition_{str(acq_num)}/camera_{camera_SN}/'
    #general_folder = f'G:/Grenoble/{date}/manip_relation_dispersion/Acquisition_{str(acq_num)}/camera_{camera_SN}/'

    folder_video_demod = general_folder + str(f_exc) +'Hz_'+ str(freq_acq) +'Hz/matData/video_demod_W'+ str(W) +'_Dt' +str(Dt)
    matfile_path = folder_video_demod + '/figdata_complex.mat'
    print('data loading...')
    Data_demod = load_complex_field(matfile_path)
    print('data loaded!')
    complex_field = np.copy(Data_demod.data)
    if plot:
        plt.figure()
        plt.imshow(np.real(complex_field)/np.max(np.abs(complex_field)))
        plt.show()
    ## process
    nb_periods = 4
    nb_pts_per_period = 256
    #index_profile_line = 34 # pour piv du 25/11 c'est 34 (W=32), pour W=64 on met donc 17
    m = nb_periods*nb_pts_per_period
    matrix = np.zeros((m,complex_field.shape[1]))
    for i in range(m):
        img = np.real(complex_field*np.exp(-2*np.pi*1j*i/nb_pts_per_period))#/np.max(np.abs(complex_field))
        img = np.real(complex_field*np.exp(-2*np.pi*1j*i/nb_pts_per_period))/np.abs(complex_field)
        if dosavgol==True:
            img = savgol_filter(img,21,3)
        #img = np.real(complex_field*np.exp(2*np.pi*1j*i/nb_pts_per_period))/np.abs(complex_field)
        p = img[index_profile_line,:]
        #plt.imshow(img)
        #plt.pause(0.01)
        matrix[i,:] = p


    # calcul du spatio-temporel du champ démodulé en temps


    X_MAXS = np.arange(matrix.shape[1])
    Y_MAXS = []
    
    if direction==-1:
        matrix = np.flip(matrix,axis=0)
    elif direction==1:
        pass
    print('shape matrix : ',matrix.shape)
    for i in range(matrix.shape[1]):
        indices = find_peaks(matrix[:,i])[0]
        i0 = indices[0]
        i1 = indices[1]
        if len(Y_MAXS)==0:
            Y_MAXS.append(i0)
        #elif ((len(Y_MAXS)!=0)&(abs(Y_MAXS[-1]-i0)<=nb_pts_per_period/3)):
        #    Y_MAXS.append(i0)
        #elif ((len(Y_MAXS)!=0)&(abs(Y_MAXS[-1]-i0)>nb_pts_per_period/3)):
        #    Y_MAXS.append(i1)
        
        elif ((len(Y_MAXS)!=0)&(abs(Y_MAXS[-1]-i0)<=abs(Y_MAXS[-1]-i1))):
            Y_MAXS.append(i0)
        elif ((len(Y_MAXS)!=0)&(abs(Y_MAXS[-1]-i0)>abs(Y_MAXS[-1]-i1))):
            Y_MAXS.append(i1)
        else:
            Y_MAXS.append(np.nan)

    #Y_MAXS = np.unwrap(Y_MAXS,discont=48)    
    ## plot the line to fit on top of the matrix imshow (test)
    if plot:
        plt.figure()
        plt.title('$f_{exc}=$'+str(f_exc)+' Hz')
        plt.imshow(matrix.T,aspect='auto')
        plt.plot(Y_MAXS,np.arange(len(matrix[0,:])),'r')
        plt.colorbar()

        plt.ylabel('position along profile [PIV box units]',fontsize=15)
        plt.xlabel('time [256 px = T = 1/$f_{exc}$]',fontsize=15)
        #plt.ylim(np.nanmin(Y_MAXS),np.nanmax(Y_MAXS))
        plt.ylim(0,55)
        
        plt.show()

    ##  fit the line on the spatio-temporal diagram

    x_px = np.arange(xlim_fit[0],xlim_fit[1],1)
    x_meters = x_px * W/2 * (dcm*1e-2)/dpx

    #x_meters = np.unwrap(x_meters)

    t_px = Y_MAXS[xlim_fit[0]:xlim_fit[1]]
    print(t_px)
    t_sec = np.array(t_px)*(1/f_exc)/nb_pts_per_period


    def linear(x,a,b):
        return a*x+b
    t_fit = np.linspace(np.min(t_sec),np.max(t_sec),100)
    opt,pcov = curve_fit(linear,t_sec,x_meters)
    if plot:
        plt.plot(t_sec,x_meters,'o')
        plt.plot(t_fit,linear(t_fit,opt[0],opt[1]))
        plt.xlabel('Time (seconds)',fontsize=15)
        plt.ylabel('Position of the wave front (meters)',fontsize=15)
        plt.title('Slope given by the fit : $v_{\phi}$ = '+str(np.round(opt[0],2))+' +- '+str(np.round(np.sqrt(pcov[0][0]),2))+' m/s',fontsize=15)
        plt.show()
    print('phase velocity = '+str(opt[0])+' +- '+str(np.sqrt(pcov[0][0]))+' m/s')
    # il faut bien prendre sqrt(pcov) pour prendre l'ecart type et non la variance !!!
    return opt[0],np.sqrt(pcov[0][0])

###########################################
# fonction listallfreq

def list_all_freq(general_folder):
    ll = os.listdir(general_folder)
    ld = []
    list_f_exc = []
    list_freq_acq = []
    for l in ll:
        if ('Hz_' in l)&(l[-2:]=='Hz'):
            ld.append(l)
            print(l)
            list_f_exc.append(int(l[:l.find('Hz')]))
            if '.' in l[l.find('_')+1:-2]:
                list_freq_acq.append(float(l[l.find('_')+1:-2]))
            elif ('.' in l[l.find('_')+1:-2])==False:
                list_freq_acq.append(int(l[l.find('_')+1:-2]))

    tab_f_exc = np.array(list_f_exc)
    tab_freq_acq = np.array(list_freq_acq,dtype='object')
    indsort = np.argsort(tab_f_exc)
    tab_f_exc = tab_f_exc[indsort]
    tab_freq_acq = tab_freq_acq[indsort]

    return tab_f_exc,tab_freq_acq