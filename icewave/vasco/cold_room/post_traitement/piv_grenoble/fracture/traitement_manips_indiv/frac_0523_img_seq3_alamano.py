#%%
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pickle
#%%
date = '0523'
name_frac_file = 'img_seq3'
#camera_SN = '22458101'

f_exc = 0.9597
freq_acq = 20

ypix_surf = 380




#%%

# ici on va juste voir le resultat que ca donne avec image j, en regardant l'amplitude 
# à plusieurs points pour estimer la courbure
frame_frac = 4603

frame1 = 4572
frame2 = 4582

x1_vals = np.array([1152,833,587,1289,1053,724,942,791,884,848,823,784])
x2_vals = np.array([1154,837,586,1290,1053,724,944,789,884,849,823,783])

y1_vals = np.array([283,285,291,289,288,280,295,310,312,283,306,304])
y2_vals = np.array([271,264,277,281,273,261,276,290,292,262,285,284])

x2fit = (x1_vals+x2_vals)/2
y2fit = (y1_vals-y2_vals)/2

indices2fit = np.where((x2fit>500)&(x2fit<1200))

a,b,c = np.polyfit(x2_vals[indices2fit],y2fit[indices2fit],2)


xth=np.linspace(np.min(x2fit),np.max(x2fit))
plt.figure()
plt.plot((x1_vals+x2_vals)/2, (y1_vals-y2_vals)/2,'o',label='(y1_vals-y2_vals)/2')
plt.plot(xth,a*xth**2+b*xth+c)
plt.ylim(0,11)
plt.legend()
plt.show()


dcmperdpx = 0.065
xth=np.linspace(np.min(x2fit)*dcmperdpx,np.max(x2fit)*dcmperdpx) * 1e-2

a,b,c = np.polyfit(x2_vals[indices2fit] * dcmperdpx * 1e-2,y2fit[indices2fit] * dcmperdpx * 1e-2,2)

plt.figure()
plt.plot(1e-2 * dcmperdpx * (x1_vals+x2_vals)/2,1e-2 * dcmperdpx * (y1_vals-y2_vals)/2,'o',label='(y1_vals-y2_vals)/2')
plt.plot(xth,a*xth**2+b*xth+c)
#plt.ylim(0,11*1e-2)
lambda_approx = 1.4
popt,pcov = curve_fit((lambda x,A,phi:A*np.cos((2*np.pi/lambda_approx)*x+phi)),x2_vals * dcmperdpx * 1e-2,y2fit * dcmperdpx * 1e-2)
plt.plot(xth,popt[0]*np.cos((2*np.pi/lambda_approx)*xth+popt[1]))
plt.legend()
plt.show()



kappa = 2*a

kappa_c_vals = np.array([kappa,popt[0]*(2*np.pi/lambda_approx)**2])

# correction angulaire (ça change un peu)
kappa_c_vals = kappa_c_vals/np.cos(np.radians(8))

kappa_c_avg = np.mean(kappa_c_vals)
kappa_c_std = np.std(kappa_c_vals)

#print(kappa)
#print(popt[0]*(2*np.pi/lambda_approx)**2)
print(kappa_c_vals)
print('kappac=',kappa_c_avg,' +- ',kappa_c_std,' m^-1')



# %% enregistrement des données
file_dict_results = 'R:/Gre25/Summary/fracture_postprocessing/resultats/fracture_results.pkl'

if os.path.exists(file_dict_results):
    with open(file_dict_results, 'rb') as f:
        dict_results = pickle.load(f)
else:
    dict_results = {}



if date in dict_results:
    dict_results[date][name_frac_file] = {'method':'ImageJ','kappa_c_vals':kappa_c_vals, 'ypix':285, 'ypix_surf':ypix_surf,'f_exc':f_exc,'frame1':frame1,'frame2':frame2}
else:
    dict_results[date] = {name_frac_file: {'method':'ImageJ', 'kappa_c_vals':kappa_c_vals, 'ypix':285, 'ypix_surf':ypix_surf,'f_exc':f_exc,'frame1':frame1,'frame2':frame2}}






ee = input('Are you sure you want to erase last results ?(y/n)')
if ee=='y':
    with open(file_dict_results, 'wb') as handle:
        pickle.dump(dict_results, handle, protocol=pickle.HIGHEST_PROTOCOL)
else:
    pass

# %%
