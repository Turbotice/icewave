# -*- coding: utf-8 -*-
"""
Created on Fri Apr  8 12:31:10 2022

@author: Turbots
"""

import cv2 
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import convolve2d
from skimage import filters
import os
from scipy.optimize import curve_fit
from scipy import stats
from PIL import Image
import pickle
from skimage.morphology import skeletonize


'''
Fonctions de démarrage
'''

def import_images(loc, nom_exp, exp_type, nom_fich = "\image_sequence\\"):
    fichiers = []                             
    liste_images = []

    fichiers = os.listdir(loc)        

    for j in range (len (fichiers)):
        if nom_exp == fichiers[j][8:13] :
            if exp_type in fichiers[j]:
                titre_exp = fichiers[j]
                path_images = str(loc + fichiers [j] + nom_fich)

    liste_images = os.listdir(path_images)

    print (path_images)
      
    #Créé un répertoire pour les résultats
         
    # if os.path.isdir(path_images[:-15] + "resultats") == False:
    #       os.mkdir(path_images[:-15] + "resultats")
          
    return path_images, liste_images,titre_exp

def import_calibration(titre_exp, date):
    if "LAS" in titre_exp:
        if float(date) >= 220412:
            mmparpixelx = 0.0974506899508849 #profondeur
            mmparpixely = 0.0390800554936788 #vertical
            mmparpixelz = 0.0421864387474003 #horizontal
            mmparpixel = 1
            if float(date) >= 220512:
                mmparpixelx = 0.150384343422359 #profondeur
                mmparpixely = 0.043816976811791 #vertical
                mmparpixelz = 0.045110366648328 #horizontal
                mmparpixel = 1
                if float(date) >= 220523:
                    mmparpixelx = 0.150811071469620 #profondeur
                    mmparpixely = 0.042888820001944 #vertical
                    mmparpixelz = 0.044245121848119 #horizontal
                    mmparpixel = 1
                    if float(date) >= 220607:
                        mmparpixelz = 0.045008078623617
                        mmparpixely = 0.047095115932238
                        mmparpixel = 1
                        if float(date) >= 221006:
                            mmparpixel = 0.28
                            mmparpixelz = 0.28
                            mmparpixely = 0.28
                            if float(date) >= 221011:
                                mmparpixel = 0.29003161344586
                                mmparpixelz = 0.29003161344586
                                mmparpixely = 0.29003161344586
                                if float(date) >= 221012:
                                    mmparpixel = 0.2001
                                    mmparpixelz = 0.2001
                                    mmparpixely = 0.2001
                                    if float(date) >= 221024:
                                        mmparpixel = 0.19434
                                        mmparpixelz = 0.19434
                                        mmparpixely = 0.19434
                                        if float(date) >= 221116:
                                            mmparpixel = 0.18763
                                            mmparpixelz = 0.18763
                                            mmparpixely = 0.18763
                                            if float(date) >= 221128:
                                                mmparpixel = 0.1803036313
                                                mmparpixelz = 0.1803036313
                                                mmparpixely = 0.1803036313
                                                if float(date) >= 221214:
                                                    mmparpixel = 0.1823586
                                                    mmparpixelz = 0.1823586
                                                    mmparpixely = 0.1823586
                                                    if float(date) >= 230103:
                                                        mmparpixel = 0.43071886979
                                                        mmparpixelz = 0.43071886979
                                                        mmparpixely = 0.43071886979
                                                        if float(date) >= 230105:
                                                            mmparpixel = 0.42824718
                                                            mmparpixelz = 0.42824718
                                                            mmparpixely = 0.42824718
                                                            if float(date) >= 230110:
                                                                mmparpixel = 0.431276146
                                                                mmparpixelz = 0.431276146
                                                                mmparpixely = 0.431276146
                                                                if float(date) >= 230117:
                                                                    mmparpixel = 0.429258241758
                                                                    mmparpixelz = 0.429258241758
                                                                    mmparpixely = 0.429258241758
                                                                    if float(date) >= 230220:
                                                                        mmparpixel = 0.41809515845806505
                                                                        mmparpixelz = 0.41809515845806505
                                                                        mmparpixely = 0.41809515845806505
                                                            
                                                            
                                                            
                                                        
                                                       
                                                    
                                                
                                               
                                        
                                


                        
        angle_cam_LAS = np.arccos(mmparpixely/mmparpixelz) * 180 / np.pi
        return mmparpixelx, mmparpixely, mmparpixelz, angle_cam_LAS, mmparpixel

    if "FSD" in titre_exp or "PIV" in titre_exp :

        if float(date) >= 220405:
            mmparpixel = 0.04                   #220405 f = 50mm
            if float(date) >= 220407:
                mmparpixel = 0.2196122526067974     #220407 f = 12mm salle 235
                if float(date) >= 220512:
                    mmparpixel = 0.230629922620838      #220512 f = 12mm salle 223
                    if float(date) >= 220516:
                        mmparpixel = 0.230578001345791      #220516 f = 12mm salle 223
                        if float(date) >= 220523:
                            mmparpixel = 0.230433333578193      #220523 f = 12mm salle 223
                            if float(date) >= 220607:
                                mmparpixel = 0.230340502325192
                                if float(date) >= 221006:
                                    mmparpixel = 0.28
                                    if float(date) >= 221011:
                                        mmparpixel = 0.29003161344586
                                        if float(date) >= 221012:
                                            mmparpixel = 0.2001
                                            if float(date) >= 221024:
                                                mmparpixel = 0.19434
                                                if float(date) >= 221116:
                                                    mmparpixel = 0.18763
                                                    if float(date) >= 221128:
                                                        mmparpixel = 0.1803036313
                                                        if float(date) >= 221214:
                                                            mmparpixel = 0.1823586
                                                            if float(date) >= 230103:
                                                                mmparpixel = 0.43071886979
                                                                if float(date) >= 230105:
                                                                    mmparpixel = 0.42824718
                                                                    if float(date) >= 230110:
                                                                        mmparpixel = 0.431276146
                                                                        if float(date) >= 230110:
                                                                            mmparpixel = 0.431276146
                                                                            if float(date) >= 230117:
                                                                                mmparpixel = 0.429258241758
                                                                                if float(date) >= 230220:
                                                                                    mmparpixel =  0.41809515845806505
                                                                                
                                                                                
                                                                                
                                                                               
                                                                        
                                                                        
                                                                       
                                                            
                                                            
                                                    
        

        return mmparpixel
    
    if  "IND" in titre_exp:
        if float(date) >= 220701 : #TIPP1
             mmparpixely = (0.03305009402751750828731107740002 + 0.03298185668063997994703113817089) /2 #horizontal
             mmparpixelz = (4.175626168548983268432967670966E-2 + 0.0416406412658754944826150322715) / 2  #vertical
             if float(date) >= 220708:
                 if "IJSP2" in titre_exp :
                     mmparpixely = 0.0316139
                     mmparpixelz = (0.0428535 + 0.04309175 + 0.0432496) / 3
                 else :
                     mmparpixely = 0.031963
                     mmparpixels = 0.042547
                     
                    
            
        
        
        
        
        
        return mmparpixely, mmparpixelz
             
    
'''
Fonctions pour craks
'''


def weighted_avg_and_std(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    average = np.average(values, weights=weights)
    variance = np.average((values-average)**2, weights=weights)
    return (average, np.sqrt(variance))


def RemoveNodes(img):
    """
    Retire les pixels qui ont plus de 3 voisins d'une image skeletonisée 
    ie Retire les noeuds des branches des fissures
    Retourne la même image sans les noeuds
    """
    output = np.copy(img)
    img2 = img/255
    
    kernel1 = np.array([[1, 1, 1],
                        [1, 0, 1],
                        [1, 1, 1]])
    
    convolve = cv2.filter2D(src=img2, ddepth=-1, kernel=kernel1)
    nodes = np.where(convolve >= 3)
    output[nodes] = 0
    
    return(output)


def cracks(img, crack_lenght_min):
    """
    Prend en entré une image skeletonisée
    et la longueur minimale en pixel des fissures que l'on prend en compte
    
    Retourne:
    liste_arg: liste où est listé chaque fissure avec ses coordonnées cartésiennes (en pixel)
    use: liste des fissuures avec longueur > crack_lenght_min
    NX, NY: taille de l'image
    """

    img = RemoveNodes(img)

    cc = cv2.connectedComponents(img,connectivity=8)


    liste_arg = [[] for i in range(0,cc[0])]
    NY,NX = cc[1].shape

    for j in range(NX):
        for i in range(NY):
            liste_arg[cc[1][i,j]].append([i,j])

    del liste_arg[0]
 
    list_big_cracks = []
    for i,list_coord in enumerate(liste_arg):
        if len(list_coord) >= crack_lenght_min:
            list_big_cracks.append(list_coord)
    
    return liste_arg, list_big_cracks, NX, NY
      
        
def analyze_cracks(liste_arg, list_big_cracks, original = None, display = False):
    """
    Input:
    liste_arg: liste où est listé chaque fissure avec ses coordonnées cartésiennes (en pixel)
    list_big_cracks: liste des fissures que l'on va considérer
    plot: if True: trace les droites associées aux fit des fissures (optionnel)
    original: chemin de l'image original pour pouvoir supperposer les droites par dessus (optionnel)
    
    Output:
    angles: array contenant l'angle d'orientation en degré de chaque fissure
    weights: array contenant la longueur en pixel de chaque fissure
    
    """
    Ncomp = len(list_big_cracks)

    angles = np.zeros(Ncomp)
    weights = np.zeros(Ncomp)
    fissures = []
    for i in range (len(list_big_cracks)):
        y = np.array (list_big_cracks[i])[:,1]
        x = np.array (list_big_cracks[i])[:,0]
        weights[i] = len(x)
        try:
            [a,b] = np.polyfit(x,y, 1)
            angles[i] = np.arctan(a)*180/np.pi
            fissures.append([a,b])
                
        except:  #les lignes parfaitement horizontales renvoient une erreur (pente infinie)
            angles[i] = 90.
    
    return angles, weights, fissures


def histo_angles(angles, weights, color=None):
    """
    Trace l'histogramme des angles pondéré par la taille des fissures
    """
    histo = plt.hist(angles,weights=weights,density=True,bins = 50, range=(-90,90),color=color)
    maxi = np.max(histo[0])
    plt.plot([-45,-45],[0,maxi],color='black',linestyle='dotted',alpha=0.7)
    plt.plot([45,45],[0,maxi],color='black',linestyle='dotted',alpha=0.7)
    plt.plot([0,0],[0,maxi],color='black',linestyle='dashed',alpha=0.7)
    plt.xlabel('angle (°)')
    plt.ylabel('densité')


def gaus(x,a1,a2,x1,x2,sigma1,sigma2):
  
  # Fit les angles avec 2 gaussiennes
  # y = histo[0]
  # x = (histo[1][:-1]+ histo[1][1:])/2
 
  # neg = angles[np.where(angles<0)]
  # negweights = weights[np.where(angles<0)]
 
  # pos = angles[np.where(angles>0)]
  # posweights = weights[np.where(angles>0)]
 
  # s1 = weighted_avg_and_std(neg,negweights)[1]
  # s2 = weighted_avg_and_std(pos,weights=posweights)[1]
    
  # popt,pcov = curve_fit(gaus,x,y,p0=[0.5,0.5,np.average(neg,weights=negweights),np.average(pos,weights=posweights),s1,s2],bounds=((0, 0, -90, 0, 0, 0), (np.inf, np.inf, 0, 90, np.inf, np.inf)))
 
  # x1 = popt[2]
  # x2 = popt[3]
  # s1 = popt[4]
  # s2 = popt[5]
 
  # x = np.linspace(-90,90,180)
  # plt.plot(x,gaus(x,*popt),color='black')
 
  # plt.title(r'$\theta_-=$'+str(int(np.round(x1,0)))+'$\pm$'+str(int(np.round(s1,0)))+'°      '+r'$\theta_+=$'+str(int(np.round(x2,0)))+'$\pm$'+str(int(np.round(s2,0)))+'°')

    return a1*np.exp(-(x-x1)**2/(2*sigma1**2))+a2*np.exp(-(x-x2)**2/(2*sigma2**2))

'''
Fonctions pour FSD BAPT
'''

# Fonction qui erode puis dilate l'image pour la rendre plus propre et faciliter la détection de contours

def erodedilate(image, kernel_iteration, kernel_size, save = False, path_images = '', name_save = ''):
    
    kernel = np.ones((kernel_size,kernel_size), np.uint8)
    
    img_erosion = cv2.erode(image, kernel, iterations=kernel_iteration)
        
    img_dilatation = cv2.dilate(img_erosion, kernel, iterations=kernel_iteration)
    
    if save :
        plt.imsave(path_images[:-15] + "resultats" + "/" + name_save + "_erodilate.tiff" , img_dilatation, cmap=plt.cm.gray )
        
    return img_dilatation
    

def conv2(x, y, mode='same'):
    return np.rot90(convolve2d(np.rot90(x, 2), np.rot90(y, 2), mode=mode), 2)


'''
Fonctions pour LAS
'''

def import_angle (date, nom_exp, loc, display = False):
    grossissement = 7.316
    path, liste, titre = import_images(loc,nom_exp,"LAS",nom_fich = '\\references')

    dx_dh = np.loadtxt(path + '\\angle_laser.txt')

    dx = dx_dh[:,0]
    dh = dx_dh[:,1]
    
    if display :
        plt.figure()
        plt.plot(dx,dh,'ko')
        

    hprime = np.polyfit(dx,dh,1)

    dh_prime = dh - hprime[1]
    xxx = np.linspace(0,max(dx),200)
    
    if display :
        plt.plot(xxx, hprime[0] * xxx + hprime[1])
        plt.xlabel("Distance laser")
        plt.ylabel("h mesuré")
    
    tab_angles = np.arctan(dh_prime/dx) * 180 / np.pi

    angle = np.mean(tab_angles)
    xxx = np.linspace(1,len(dh_prime) + 1, len(dh_prime))
    er_angle = (max(tab_angles) - min(tab_angles))/2
    if display :
        plt.figure()
        plt.plot(dx, dh_prime, 'mx')
        plt.figure()
        plt.plot(xxx, tab_angles, 'ko')
        plt.xlabel("Mesure")
        plt.ylabel("Angle (en radian)")
        plt.title('Angles')
        print('angle moyen : ', angle, 'erreur', er_angle )
        plt.figure()
        plt.plot(dx, tab_angles, 'ko')
        plt.xlabel("distance laser")
        plt.ylabel("Angle (en radian)")
    grossissement = 1 / np.tan(angle* np.pi / 180)
    grossissements = 1 / np.tan(tab_angles* np.pi / 180)
    er_grossissement = (max(grossissements) - min (grossissements))/2
    if display :
        print('grossissement :',grossissement, "erreur ", er_grossissement)
    
    return grossissement, er_grossissement, angle, er_angle


'''
Fonctions pour FFT
'''

def find_best_lin (Y, X = False, range_taille = 3, pas_taille = 12):
    #prend une courbe et trouve la meilleure zone avec un fit linéaire, zone de taille 1 / 1 + range_taille et on decale de 1/pas_taille à chaque fois
     
    R = []
    length = len (Y)
    
    for taille_test in range (1 , range_taille + 1) :
        for pos_test in range (pas_taille) :
            if 1 >= (1/taille_test + pos_test/pas_taille):
                 
                zone_test = Y [int(length/ pas_taille * pos_test) : int(length *  (1/taille_test + pos_test/pas_taille) ) ]
                if np.mean(X) == False :
                    x = np.arange(len (zone_test))
                else :
                    x = X [int(length/ pas_taille * pos_test) : int(length *  (1/taille_test + pos_test/pas_taille) ) ]
                
                slope, intercept, r_value, p_value, std_error = stats.linregress(x,zone_test)
                r = r_value ** 2
                R.append ([[slope,intercept], abs(r), taille_test, pos_test])
                
    R = np.asarray(R)    
    max_r = max(R[:,1]) #coeff de corélation le meilleur
    max_corel = R[np.argmax(R[:,1])]   
    p = max_corel[0] #[pente, origine] où il y a le meilleur coeff correl
    
    best_taille = max_corel[2] #taille de la zone la meilleure
    best_pos = max_corel[3] #pos de la meilleure zone
    
    #return pente et origine, coeff correlation, pos, taille de la zone fittée
    
    return p, max_r, best_taille, best_pos
    
    
def demodulation_gld_sp_ba(t,s, fexc, t1):
    # s signal de l'elevation (x,t)
    # # Exponentielle complexe pour la demodulation:
    c = np.mean(s*np.exp(-1j * 2 * np.pi * t[None,:] * fexc),axis=1)
    etademod = np.real(c[:,None]*np.exp(1j*2*np.pi*t1[ None,:]))
    return c, t1, etademod

'''
Fonctions pour le dico
'''
def open_dico():
    
    a_file = open("D:\Banquise\Baptiste\Resultats\\Dictionnaire.pkl", "rb")
    
    dico = pickle.load(a_file)
    
    a_file.close
    
    return dico
    
def save_dico(dico):
    
    a_file = open("D:\Banquise\Baptiste\Resultats\\Dictionnaire.pkl", "wb")
    
    pickle.dump(dico, a_file)
    
    a_file.close()
    

def add_dico(dico, date, nom_exp, name, value):
    dico[date][nom_exp][name] = value
    return dico
    
    

def remove_dico(dico, date, nom_exp, name):
    del dico[date][nom_exp][name]
    
def rename_variable_dico (dico,date, nom_exp, old_name, new_name):
    value = dico[date][nom_exp][old_name]
    dico = add_dico(dico, date, nom_exp, new_name, value)
    remove_dico(dico,date,nom_exp,old_name)
    return dico
    
        
    