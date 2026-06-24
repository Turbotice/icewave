#%%
import numpy as np
import matplotlib.pyplot as plt
import pickle
import os
import sys

icewave_path = 'C:/Users/Vasco Zanchi/Documents/git_turbotice/icewave/icewave/'
sys.path.append(icewave_path)

import tools.rw_data as rw

from vasco.tools.clickonfigures import profile_line_on_image_2clicks
from vasco.tools.clickonfigures import get_n_points_onimage

from functions_fracture_analysis import click_on_fracture_path_plot_time_evol, click2extract_amplitude, plot_elevation_refnotbroken_and_broken

#%%
#disk = 'L:'# disk is Elements on adour
disk = 'C:'# disk is Elements on adour

date = '0211'

#data_path = f'{disk}/Share_hublot/Data'
data_path = f'C:/Users/Vasco Zanchi/Desktop/Saguenay2024'
daily_drone_data_path = f'{data_path}/{date}/Drones'

velocity_field_path = f'{daily_drone_data_path}/exact_solution_real_field_stereo_0211_2024_rectangular_grid.h5'

# traitement stéphane : postitions des fractures (obtenue avec divergence du champ de vitesses)
fractures_positions_path = f'{daily_drone_data_path}/Results/fracture_positions.pkl'

#%% chargement des positions des fractures détectées
with open(fractures_positions_path, 'rb') as file:
    fractures_positions_data = pickle.load(file)


%matplotlib inline

arr_bools = np.bool_(fractures_positions_data['binary']==1)
data2show = fractures_positions_data['tmaxs_kappa2_t'].astype(float).copy()
data2show[~arr_bools] = np.nan
plt.figure()
plt.imshow(data2show)
plt.show()





# fabriquer le jeu de donnée d'entrée pour clustering

yindices_vals = []
xindices_vals = []
tfrac_sec_vals = []

for i in range(data2show.shape[0]):
    for j in range(data2show.shape[1]):
        if not(np.isnan(data2show[i,j])):
            yindices_vals.append(i)
            xindices_vals.append(j)
            tfrac_sec_vals.append(data2show[i,j])

yindices_vals = np.array(yindices_vals)
xindices_vals = np.array(xindices_vals)
tfrac_sec_vals = np.array(tfrac_sec_vals)

X = np.zeros((len(yindices_vals), 3))
X[:,0] = yindices_vals
X[:,1] = xindices_vals
X[:,2] = tfrac_sec_vals







#%%
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.datasets import make_blobs
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import SpectralClustering
from sklearn.cluster import DBSCAN
import hdbscan
from sklearn.cluster import AgglomerativeClustering


# ── 1. Generate sample data ──────────────────────────────────────────────────
#X, y_true = make_blobs(n_samples=300, centers=4, cluster_std=0.8, random_state=42)
X = StandardScaler().fit_transform(X)

# ── 2. Fit K-Means ───────────────────────────────────────────────────────────
n_clusters = 100


sc = SpectralClustering(n_clusters=n_clusters, affinity="nearest_neighbors")
#db = DBSCAN(eps=0.5, min_samples=4)
#clusterer = hdbscan.HDBSCAN(min_cluster_size=10)
labels = sc.fit_predict(X)

#model = AgglomerativeClustering(n_clusters=50, linkage='single')
#labels = model.fit_predict(X)

k = 100
#kmeans = KMeans(n_clusters=k, init="k-means++", n_init=10, random_state=42)
#kmeans.fit(X)

#labels     = kmeans.labels_
#centers    = kmeans.cluster_centers_
#inertia    = kmeans.inertia_





"""print(f"Inertia (WCSS): {inertia:.2f}")
# ── 3. Elbow Method – find optimal k ─────────────────────────────────────────
inertias = []
K_range  = range(1, 11)

for k in K_range:
    km = KMeans(n_clusters=k, init="k-means++", n_init=10, random_state=42)
    km.fit(X)
    inertias.append(km.inertia_)
"""
# ── 4. Plot ───────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 1, figsize=(10, 7))

# Cluster plot
colors = plt.cm.tab10(np.linspace(0, 1, 4))
for i in range(k):

    mask = labels == i
    axes.scatter(X[mask, 0], X[mask, 1], s=40, label=f"Cluster {i+1}")
#axes[0].scatter(centers[:, 0], centers[:, 1],
#                c="black", s=200, marker="X", zorder=5, label="Centroids")
#axes.set_title(f"K-Means Clustering (k={k})")
axes.set_xlabel("Feature 1")
axes.set_ylabel("Feature 2")
#axes[0].legend()

"""# Elbow curve
axes[1].plot(K_range, inertias, "bo-", linewidth=2, markersize=8)
axes[1].axvline(x=4, color="red", linestyle="--", label="Optimal k=4")
axes[1].set_title("Elbow Method")
axes[1].set_xlabel("Number of Clusters (k)")
axes[1].set_ylabel("Inertia (WCSS)")
axes[1].legend()"""

plt.tight_layout()
#plt.savefig("kmeans_result.png", dpi=150)
plt.show()


