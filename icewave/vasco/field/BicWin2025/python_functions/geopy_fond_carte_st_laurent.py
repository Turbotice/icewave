#%%
%matplotlib qt
import geopandas as gpd
import matplotlib.pyplot as plt
import contextily as ctx
from shapely.geometry import Point

# Coordonnées du centre de baie du ha!ha!
center_lon, center_lat = -68.818, 48.348

# Exemple de points
points = [
    (-68.818, 43.348)
]

# Création du GeoDataFrame
gdf = gpd.GeoDataFrame(
    geometry=[Point(lon, lat) for lon, lat in points],
    crs="EPSG:4326"
).to_crs(epsg=3857)

# Figure haute résolution
fig, ax = plt.subplots(figsize=(10, 10), dpi=200)

# Tracer les points
gdf.plot(ax=ax, color='red', markersize=50)

# Définir la zone d’affichage (~10 km)
center = gpd.GeoSeries([Point(center_lon, center_lat)], crs="EPSG:4326").to_crs(epsg=3857).iloc[0]
delta = 5000  # en mètres
ax.set_xlim(center.x - delta, center.x + delta)
ax.set_ylim(center.y - delta, center.y + delta)

# Ajouter le fond de carte avec un zoom plus élevé
ctx.add_basemap(
    ax,
    source=ctx.providers.OpenStreetMap.Mapnik,
    zoom=13,  # augmentez à 14–15 pour plus de détails
    alpha=1
)

ax.set_title("Baie du Ha!Ha! (haute résolution)", fontsize=14)
plt.tight_layout()
plt.show()
µ# %%

# %%
